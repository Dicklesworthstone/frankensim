//! The original molecular/circular-wall law in the continuous modal field.
//! Geometry, quadrature, radial projection and wave-chart propagation retain
//! their original owners. Molecular Cp and wall heat exchange share the same
//! constitutive response as MolecularDuct; their attenuation is never added.

use super::*;
use crate::molecular::{MolecularDuctError, MolecularWallResolution, gas_key, local, state_key};
use fs_material::relaxation::{RelaxationAbsorption, RelaxationState, RelaxingMedium};
use fs_phs::{ZkDiffusive, ZkDiffusiveBasis};
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
struct MaterialPoint {
    radius: f64,
    x_m: f64,
    state: usize,
    wall: Option<ZkDiffusive>,
}

#[derive(Clone, Debug)]
struct MolecularStep {
    h: f64,
    left_m: f64,
    positions_m: [f64; 2],
    points: [usize; 2],
    flares: [f64; 2],
}

#[derive(Clone, Debug)]
struct Binding {
    original: ModalGasFieldProfile,
    law: RelaxationAbsorption,
    wall: Option<MolecularWallResolution>,
    states: Vec<RelaxationState>,
    points: Vec<MaterialPoint>,
    steps: Vec<MolecularStep>,
    outlet_point: usize,
}

/// Source-bound molecular dispersion in every original retained radial mode.
///
/// This is the existing axisymmetric, frozen-background modal approximation.
/// Without wall diffusion, original molecular bulk compliance is used directly.
/// With a wall, the native coupled circular heat/shear response supplies the
/// original plane-ZK higher-mode closure; this is not a resolved higher-mode
/// thermoviscous boundary layer or an evolving heat/species simulation.
/// Frequency preparation shares immutable material/geometry storage. Repeated
/// physical source drives use the existing prepared field without re-evaluation.
#[derive(Clone, Debug)]
pub struct MolecularModalGasFieldProfile {
    binding: Arc<Binding>,
}

impl ModalGasFieldProfile {
    /// Bind the original relaxation law at the exact cached physical gas-field
    /// quadrature and outlet. None selects adiabatic slip walls; Some binds
    /// the original positive circular diffusion measure to every gas/radius.
    /// Its pole count, admitted band and sampled error are independent of the
    /// original axial steps and retained radial-mode count.
    ///
    /// Exact state/radius caches share frequency evaluations and one dimensionless
    /// diffusion eigensystem, without rounding coordinates or merging processes.
    /// An empty law without walls delegates the original lossless owner exactly.
    /// Source, diffusion, geometry and cancellation refusals publish no binding.
    pub fn with_molecular_relaxation(
        self,
        law: RelaxationAbsorption,
        wall: Option<MolecularWallResolution>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularModalGasFieldProfile, MolecularDuctError> {
        poll(cancelled, 0)?;
        let evaluate = |x| -> Result<RelaxationState, MolecularDuctError> {
            Ok(law.evaluate(&RelaxingMedium::from_thermochemical(&self.field.state_at(x)?)?)?)
        };
        let inlet = evaluate(0.0)?;
        let mut state_cache = BTreeMap::from([(state_key(&inlet), 0)]);
        let mut states = vec![inlet];
        let basis = wall
            .map(|r| {
                ZkDiffusiveBasis::new(r.poles_per_field).map_err(|reason| {
                    DuctError::WallDiffusionUnresolved {
                        reason,
                        radius_m: self.profile.inlet_radius(),
                        x_m: 0.0,
                    }
                })
            })
            .transpose()?;
        poll(cancelled, 0)?;
        let mut point_cache = BTreeMap::new();
        let mut wall_cache = BTreeMap::new();
        let mut points = Vec::new();
        let count = self.steps_used();
        let mut steps = Vec::with_capacity(count);
        let uniform_coordinates = uniform_step_coordinates(self.profile.length(), self.steps_per_span);
        let outlet_point = {
            let mut sample = |x_m: f64, radius: f64| -> Result<usize, MolecularDuctError> {
                let state = if self.quadrature.is_empty() {
                    0
                } else {
                    let original = evaluate(x_m)?;
                    let key = state_key(&original);
                    if let Some(&index) = state_cache.get(&key) {
                        index
                    } else {
                        let index = states.len();
                        states.push(original);
                        state_cache.insert(key, index);
                        index
                    }
                };
                let key = (state, radius.to_bits());
                if let Some(&index) = point_cache.get(&key) {
                    return Ok(index);
                }
                let gas = states[state].gas();
                let diffusion = if let Some(r) = wall {
                    let key = (radius.to_bits(), gas_key(gas));
                    if let Some(model) = wall_cache.get(&key) {
                        Some(ZkDiffusive::clone(model))
                    } else {
                        let model = basis
                            .as_ref()
                            .expect("admitted wall basis")
                            .bind(
                                radius,
                                gas.density,
                                gas.dynamic_viscosity,
                                gas.gamma,
                                gas.prandtl,
                                r.band_rad_s,
                                r.maximum_relative_error,
                            )
                            .map_err(|reason| DuctError::WallDiffusionUnresolved {
                                reason,
                                radius_m: radius,
                                x_m,
                            })?;
                        wall_cache.insert(key, model.clone());
                        Some(model)
                    }
                } else {
                    None
                };
                let index = points.len();
                points.push(MaterialPoint {
                    radius,
                    x_m,
                    state,
                    wall: diffusion,
                });
                point_cache.insert(key, index);
                Ok(index)
            };
            for j in 0..count {
                poll(cancelled, j)?;
                let (h, left_m, positions_m, radii, flares) = if self.quadrature.is_empty() {
                    let (h, left, x) = uniform_coordinates(j);
                    let radii = [self.profile.radius(x[0])?, self.profile.radius(x[1])?];
                    let flares = [
                        self.profile.slope(x[0])? / radii[0],
                        self.profile.slope(x[1])? / radii[1],
                    ];
                    (h, left, x, radii, flares)
                } else {
                    let step = &self.quadrature[j];
                    (
                        step.h,
                        step.left_m,
                        step.points.each_ref().map(|p| p.x_m),
                        step.points.each_ref().map(|p| p.radius),
                        step.points.each_ref().map(|p| p.flare),
                    )
                };
                if flares.iter().any(|x| !x.is_finite()) {
                    return Err(bad("molecular modal geometric coupling is not finite").into());
                }
                let point_indices = [
                    sample(positions_m[0], radii[0])?,
                    sample(positions_m[1], radii[1])?,
                ];
                steps.push(MolecularStep {
                    h,
                    left_m,
                    positions_m,
                    points: point_indices,
                    flares,
                });
            }
            poll(cancelled, count)?;
            sample(self.profile.length(), self.profile.outlet_radius())?
        };
        poll(cancelled, count)?;
        Ok(MolecularModalGasFieldProfile {
            binding: Arc::new(Binding {
                original: self,
                law,
                wall,
                states,
                points,
                steps,
                outlet_point,
            }),
        })
    }
}

impl MolecularModalGasFieldProfile {
    /// Original continuous geometry, primitive field and axial resolution.
    pub fn original(&self) -> &ModalGasFieldProfile {
        &self.binding.original
    }

    /// The original declared molecular law, including every species/process.
    pub fn relaxation_law(&self) -> &RelaxationAbsorption {
        &self.binding.law
    }

    /// Explicit original circular wall resolution, when requested.
    pub fn wall_resolution(&self) -> Option<MolecularWallResolution> {
        self.binding.wall
    }

    /// Exact distinct material states; index zero is the original inlet.
    pub fn material_states(&self) -> &[RelaxationState] {
        &self.binding.states
    }

    /// Original molecular material at a physical Magnus node (mouth to inlet).
    pub fn relaxation_at_quadrature(&self, step: usize, node: usize) -> Option<&RelaxationState> {
        let point = *self.binding.steps.get(step)?.points.get(node)?;
        Some(&self.binding.states[self.binding.points[point].state])
    }

    /// Original quadrature abscissa in metres, with the same ordering as above.
    pub fn quadrature_position_m(&self, step: usize, node: usize) -> Option<f64> {
        self.binding.steps.get(step)?.positions_m.get(node).copied()
    }

    /// The original outlet material used by every matched terminal mode.
    pub fn outlet_relaxation(&self) -> &RelaxationState {
        &self.binding.states[self.binding.points[self.binding.outlet_point].state]
    }

    /// Prepare the native coherent field with the same local molecular law
    /// through the last node and the complete matched outlet. Analytic radiation
    /// retains the original nondispersive outlet-gas exterior assumption;
    /// it does not silently extend molecular absorption into free space.
    pub fn prepare(
        &self,
        omega: f64,
        load: MouthLoad<'_>,
        modes: ModeCount,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedModalProfile, MolecularDuctError> {
        self.prepare_load(omega, load.into(), modes, cancelled)
    }

    /// The full original radiation matrix in an explicitly supplied exterior.
    pub fn prepare_with_flanged_mouth(
        &self,
        exterior: &GasState,
        omega: f64,
        modes: ModeCount,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedModalProfile, MolecularDuctError> {
        self.prepare_load(omega, ProfileLoad::Flanged(exterior), modes, cancelled)
    }

    /// The original complete physical mouth matrix, without diagonal reduction.
    pub fn prepare_with_mouth_matrix(
        &self,
        omega: f64,
        load: &ModalMouthImpedance,
        modes: ModeCount,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedModalProfile, MolecularDuctError> {
        poll(cancelled, 0)?;
        load.admit(
            &self.original().outlet,
            omega,
            self.original().profile.outlet_radius(),
            modes.0,
        )?;
        self.prepare_load(omega, ProfileLoad::Matrix(load), modes, cancelled)
    }

    fn prepare_load(
        &self,
        omega: f64,
        load: ProfileLoad<'_>,
        modes: ModeCount,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedModalProfile, MolecularDuctError> {
        poll(cancelled, 0)?;
        modes.validate()?;
        if !(omega.is_finite() && omega > 0.0) {
            return Err(bad("molecular modal frequency must be finite and positive").into());
        }
        let data = &self.binding;
        let original = &data.original;
        if data.wall.is_none() && data.states.iter().all(|s| s.processes().is_empty()) {
            return Ok(
                original
                    .prepare_load(omega, LossModel::Lossless, load, modes, cancelled)?
                    .retain_molecular_field(self.clone()),
            );
        }
        let mut caloric = Vec::with_capacity(data.states.len());
        for state in &data.states {
            poll(cancelled, 0)?;
            caloric.push(local::caloric(state, omega)?);
        }
        let mut waves = Vec::with_capacity(data.points.len());
        for point in &data.points {
            poll(cancelled, 0)?;
            let (cp, k) = caloric[point.state];
            waves.push(local::wave(
                &data.states[point.state],
                point.radius,
                point.x_m,
                point.wall.as_ref(),
                omega,
                cp,
                k,
            )?);
        }
        let run = propagate_resolved(
            &original.profile,
            &original.inlet,
            &original.outlet,
            omega,
            LossModel::Lossless,
            Some(waves[data.outlet_point]),
            load,
            modes.0,
            data.steps.len(),
            true,
            cancelled,
            &mut |j, basis, zref| {
                let step = &data.steps[j];
                let node = |i: usize| -> Result<Local, DuctError> {
                    let point = step.points[i];
                    generator_from_wave(
                        &waves[point],
                        basis,
                        zref,
                        data.points[point].radius,
                        step.flares[i],
                    )
                };
                Ok(PropagationStep {
                    h: step.h,
                    left_m: step.left_m,
                    local: [node(0)?, node(1)?],
                })
            },
        )?;
        Ok(
            PreparedModalProfile::from_propagated(
                &original.inlet,
                &original.outlet,
                Some(Arc::clone(&original.field)),
                omega,
                modes,
                run,
            )?
            .retain_molecular_field(self.clone()),
        )
    }
}
