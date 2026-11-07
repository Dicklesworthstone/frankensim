//! Accepted normal stresses of the original longitudinal Maxwell members.
//!
//! Each layer retains its own equilibrium modulus and every signed elastic arm
//! extension. On the native reference coordinates x=[axial, symmetric bending,
//! antisymmetric bending], axial strain is x0/L and the normal-stress gradient
//! is E*(3*beta*(1-2*xi)*x1+x2)/L^2. Each arm contributes the same expression
//! with its original delta_E and its actual elastic extensions, not current x.
//! Thus constrained memories can govern even at zero nodal displacement.
//! Reconstruction uses each original arm's native channel force: F/A for axial
//! stress, F*L/(I*beta) for symmetric bending and F*L/I for antisymmetric bending.
//! This avoids losing representable stress when the diagnostic arm extension
//! underflows. A nonzero canonical memory with an unrepresentable force refuses.
//!
//! Normal stress is affine along the member and through its cross-section.
//! Both ends and each original layer radius bound this declared interpolation.
//! Transverse-shear memory remains in the mechanics but is not a longitudinal
//! normal stress. This observer does not assess shear/principal-stress failure,
//! rupture, fatigue, or peaks between accepted temporal endpoints.

use super::PreparedFrameRelaxationTransient;
use crate::frame2d::relaxation::FrameStrainChannel;
use crate::frame2d::strength::FrameTensileError;
use fs_exec::Cx;
use fs_material::section::strength::{AxialTensileCapacity, LayerTensileAssessment};

/// One original hereditary layer at one physical member end.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameRelaxationTensileSample<'a> {
    pub element: usize,
    /// Original section layer, never an area-weighted replacement.
    pub layer: usize,
    /// Zero is the original start node, one the original end node.
    pub end: usize,
    pub capacity: &'a AxialTensileCapacity,
    pub mean_normal_stress_pa: f64,
    /// Signed local-y gradient; includes equilibrium and every original E arm.
    pub normal_stress_gradient_pa_per_m: f64,
    pub assessment: LayerTensileAssessment,
}

/// Immutable observations ordered by selected member, layer, then start/end.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameRelaxationTensileAssessment<'a> {
    pub time_s: f64,
    samples: Vec<FrameRelaxationTensileSample<'a>>,
    governing: Option<usize>,
}

impl<'a> FrameRelaxationTensileAssessment<'a> {
    pub fn samples(&self) -> &[FrameRelaxationTensileSample<'a>] {
        &self.samples
    }

    /// Largest tensile utilization; original traversal order breaks exact ties.
    pub fn governing_sample(&self) -> Option<&FrameRelaxationTensileSample<'a>> {
        self.governing.map(|index| &self.samples[index])
    }

    pub fn reaches_limit(&self) -> bool {
        self.governing_sample()
            .is_some_and(|sample| sample.assessment.tensile_utilization >= 1.0)
    }
}

impl PreparedFrameRelaxationTransient<'_> {
    /// Observe the total tensile normal stress in every selected Maxwell member.
    ///
    /// Capacities follow `source().members()` order, one exact original section
    /// per member. This covers those hereditary members, not unselected members
    /// of the containing frame. Their native preparation already refuses a
    /// simultaneous material Kelvin dashpot, so no rate stress is omitted here.
    /// Rayleigh damping is not a constitutive fibre stress.
    ///
    /// Actual accepted nodal displacement and elastic arm extensions are read
    /// together. Equilibrium and arm contributions add with their signs before
    /// extracting extrema; neither a relaxed modulus nor the net section force
    /// can substitute for these per-layer histories. Original source state,
    /// reference geometry and the Timoshenko bending fraction remain fixed.
    ///
    /// `maximum_layer_samples` bounds twice the selected layer count before any
    /// allocation. Work is linear in layers plus original memory coordinates;
    /// only the returned samples allocate. No source query, time advance, stress
    /// update or intermediate sampled field is performed. Cancellation and all
    /// refusals leave the complete native trajectory unchanged.
    pub fn assess_tensile_strength<'a>(
        &self,
        capacities: &'a [AxialTensileCapacity],
        maximum_layer_samples: usize,
        cx: &Cx<'_>,
    ) -> Result<FrameRelaxationTensileAssessment<'a>, FrameTensileError> {
        checkpoint(cx)?;
        if capacities.len() != self.source.members.len() {
            return Err(FrameTensileError::Invalid(
                "one tensile capacity per original selected Maxwell member required",
            ));
        }
        let mut count = 0_usize;
        for (member, capacity) in self.source.members.iter().zip(capacities) {
            checkpoint(cx)?;
            if member.source.elastic() != capacity.section() {
                return Err(FrameTensileError::SourceMismatch {
                    element: member.element,
                });
            }
            count = capacity
                .section()
                .layers()
                .len()
                .checked_mul(2)
                .and_then(|n| count.checked_add(n))
                .filter(|&n| n <= maximum_layer_samples)
                .ok_or(FrameTensileError::Budget)?;
        }
        let mut samples = Vec::new();
        samples
            .try_reserve_exact(count)
            .map_err(|_| FrameTensileError::Allocation)?;
        let displacement = &self.native.report().displacement;
        for (member, capacity) in self.source.members.iter().zip(capacities) {
            checkpoint(cx)?;
            let mut x = [0.0; 3];
            for (coordinate, weights) in x.iter_mut().zip(&member.weights) {
                for (&dof, &weight) in member.dofs.iter().zip(weights) {
                    *coordinate = finite(*coordinate + product(weight, displacement[dof])?)?;
                }
            }
            let length = member.length_m;
            let epsilon = quotient(x[0], length)?;
            let symmetric = product(
                3.0 * member.bending_fraction,
                quotient(quotient(x[1], length)?, length)?,
            )?;
            let antisymmetric = quotient(quotient(x[2], length)?, length)?;
            let first_sample = samples.len();
            for (layer, geometry) in capacity.section().layer_integrals().iter().enumerate() {
                checkpoint(cx)?;
                let young = geometry.young_modulus_pa();
                let mean = product(young, epsilon)?;
                for (end, sign) in [1.0, -1.0].into_iter().enumerate() {
                    samples.push(FrameRelaxationTensileSample {
                        element: member.element,
                        layer,
                        end,
                        capacity,
                        mean_normal_stress_pa: mean,
                        normal_stress_gradient_pa_per_m: product(
                            young,
                            finite(sign * symmetric + antisymmetric)?,
                        )?,
                        // Filled only after all signed arm stresses participate.
                        assessment: LayerTensileAssessment {
                            minimum_axial_stress_pa: 0.0,
                            maximum_axial_stress_pa: 0.0,
                            tensile_utilization: 0.0,
                        },
                    });
                }
            }
            for memory in &self.memory[member.coordinate_range.clone()] {
                checkpoint(cx)?;
                let coordinate = &memory.coordinate;
                if coordinate.channel == FrameStrainChannel::TransverseShear {
                    continue;
                }
                let geometry = &capacity.section().layer_integrals()[coordinate.layer];
                let force = memory.opposing_force_n;
                if force == 0.0 && self.canonical_state()[coordinate.native_state] != 0.0 {
                    return Err(unrepresentable());
                }
                let contribution = match coordinate.channel {
                    FrameStrainChannel::Axial => quotient(force, geometry.area_m2())?,
                    FrameStrainChannel::SymmetricBending => product(
                        quotient(
                            quotient(force, geometry.second_moment_m4())?,
                            member.bending_fraction,
                        )?,
                        length,
                    )?,
                    FrameStrainChannel::AntisymmetricBending => {
                        product(quotient(force, geometry.second_moment_m4())?, length)?
                    }
                    FrameStrainChannel::TransverseShear => unreachable!(),
                };
                for end in 0..2 {
                    let sample = &mut samples[first_sample + 2 * coordinate.layer + end];
                    match coordinate.channel {
                        FrameStrainChannel::Axial => {
                            sample.mean_normal_stress_pa =
                                finite(sample.mean_normal_stress_pa + contribution)?;
                        }
                        FrameStrainChannel::SymmetricBending => {
                            let signed = if end == 0 {
                                contribution
                            } else {
                                -contribution
                            };
                            sample.normal_stress_gradient_pa_per_m =
                                finite(sample.normal_stress_gradient_pa_per_m + signed)?;
                        }
                        FrameStrainChannel::AntisymmetricBending => {
                            sample.normal_stress_gradient_pa_per_m =
                                finite(sample.normal_stress_gradient_pa_per_m + contribution)?;
                        }
                        FrameStrainChannel::TransverseShear => unreachable!(),
                    }
                }
            }
        }
        let mut governing = None;
        let mut largest = -1.0;
        for (index, sample) in samples.iter_mut().enumerate() {
            checkpoint(cx)?;
            let radius = sample.capacity.section().layers()[sample.layer].outer_radius_m;
            let amplitude = product(radius, sample.normal_stress_gradient_pa_per_m.abs())?;
            sample.assessment = sample
                .capacity
                .assess_layer_normal_stress(sample.layer, sample.mean_normal_stress_pa, amplitude)
                .map_err(|source| FrameTensileError::Material {
                    element: sample.element,
                    source,
                })?;
            if sample.assessment.tensile_utilization > largest {
                largest = sample.assessment.tensile_utilization;
                governing = Some(index);
            }
        }
        checkpoint(cx)?;
        Ok(FrameRelaxationTensileAssessment {
            time_s: self.native.report().time_s,
            samples,
            governing,
        })
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), FrameTensileError> {
    cx.checkpoint().map_err(|_| FrameTensileError::Cancelled)
}

fn finite(value: f64) -> Result<f64, FrameTensileError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(unrepresentable())
    }
}

fn product(left: f64, right: f64) -> Result<f64, FrameTensileError> {
    let value = finite(left * right)?;
    if left != 0.0 && right != 0.0 && value == 0.0 {
        Err(unrepresentable())
    } else {
        Ok(value)
    }
}

fn quotient(value: f64, scale: f64) -> Result<f64, FrameTensileError> {
    let result = finite(value / scale)?;
    if value != 0.0 && result == 0.0 {
        Err(unrepresentable())
    } else {
        Ok(result)
    }
}

fn unrepresentable() -> FrameTensileError {
    FrameTensileError::Invalid("native Maxwell frame tensile stress is unrepresentable")
}
