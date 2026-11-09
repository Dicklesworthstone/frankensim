//! Original layer memories driven by a continuous coherent CST/DKT field.

use super::super::{
    C64, Cx, Dd, MaterialPlane, TensorPronyShell, TensorShellHarmonicError, checkpoint, invalid,
    section::principal_strain,
};
use super::{PreparedFacet, TensorShellMaterialBudget, cap};
use crate::PlateError;
use fs_fft::harmonic_window::{HarmonicPhaseInterval, phase_angular_interval_from_clock_dd};
use fs_material::visco::{MAX_HARMONIC_MATERIAL_LINES, TensorPronyHarmonicWindowPoint, ViscoError};
use fs_material::visco_tensor::{
    TensorPronyError, TensorPronyHistory,
    stress_control::{
        PreparedTensorMixedControl, TensorPlaneStressHarmonicBudget, TensorStressControlBudget,
    },
};
use std::sync::Arc;

const BARYCENTRIC: [[f64; 3]; 3] = [[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]];
const GAUSS: f64 = 0.5773502691896258;

/// Whole shell ceilings plus the original material owners' independent caps.
#[derive(Clone, Copy, Debug)]
pub struct TensorShellMaterialCausalBudget {
    /// Six physical quadrature points per layer. Payload also includes twelve
    /// accepted and twelve candidate histories, vertex/face gates and scratch.
    pub field: TensorShellMaterialBudget,
    pub maximum_frequencies: usize,
    pub control: TensorStressControlBudget,
    pub temporal: TensorPlaneStressHarmonicBudget,
}

struct LayerHistory {
    control: PreparedTensorMixedControl,
    virgin: TensorPronyHistory,
    faces_m: [f64; 2],
    to_principal: [[f64; 3]; 3],
}

/// Immutable original geometry and original coupled plane-stress controllers.
/// The frequency bank belongs to each interval, so a continuous drive may
/// change carriers without recreating or discarding any material memory.
pub struct TensorShellMaterialCausalObserver<'a> {
    source: &'a TensorPronyShell,
    binding: Arc<()>,
    facets: Vec<PreparedFacet>,
    layers: Vec<Vec<LayerHistory>>,
    points: usize,
    budget: TensorShellMaterialCausalBudget,
}

/// Accepted original state at an explicit normalized DD source-clock instant.
///
/// The private bank has six vertex/face and six quadrature histories per layer.
/// Starting all of them virgin and driving the original affine triangle and
/// linear thickness field preserves its bilinear spatial memory. Every interior
/// state is consequently a convex interpolation of the vertex/face states;
/// their full transient strain bounds admit the entire layer, not just samples.
/// Arbitrary history injection would break this invariant and is not exposed.
#[derive(Clone)]
pub struct TensorShellMaterialCausalHistory<'a> {
    source: &'a TensorPronyShell,
    binding: Arc<()>,
    time_s: Dd,
    full_displacement: Vec<f64>,
    histories: Vec<Vec<[TensorPronyHistory; 12]>>,
    stored_j: f64,
    work_j: f64,
    dissipated_j: f64,
}

impl<'a> TensorShellMaterialCausalHistory<'a> {
    pub const fn source(&self) -> &'a TensorPronyShell {
        self.source
    }
    pub const fn time_s_dd(&self) -> Dd {
        self.time_s
    }
    /// Full node-major translations/axial rotations, including stationary supports.
    pub fn full_displacement(&self) -> &[f64] {
        &self.full_displacement
    }
    pub const fn stored_energy_j(&self) -> f64 {
        self.stored_j
    }
    pub const fn stress_work_j(&self) -> f64 {
        self.work_j
    }
    pub const fn dissipated_energy_j(&self) -> f64 {
        self.dissipated_j
    }
}

/// Original layer quadrature receipt; endpoint and bounds use tensor shears.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorShellMaterialCausalPoint {
    pub element: usize,
    pub layer: usize,
    pub triangle_quadrature_point: usize,
    pub thickness_quadrature_point: usize,
    pub triangle_barycentric: [f64; 3],
    pub thickness_position_m: f64,
    pub volume_m3: f64,
    /// Ordered engineering [epsilon_aa, epsilon_bb, gamma_ab] in layer axes.
    pub in_plane_strain: Vec<[C64; 3]>,
    pub in_plane_bias: [f64; 3],
    /// Full original tensor-shear Voigt strain/stress, including transverse response.
    pub end_strain: [f64; 6],
    pub end_stress: [f64; 6],
    pub strain_bounds: [[f64; 6]; 2],
    pub material: TensorPronyHarmonicWindowPoint,
}

/// Independently integrated interval account on original reference volumes.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorShellMaterialCausalField {
    pub start_time_s_dd: Dd,
    pub end_time_s_dd: Dd,
    pub duration_s: f64,
    pub angular_frequencies_rad_s: Vec<f64>,
    pub points: Vec<TensorShellMaterialCausalPoint>,
    pub start_stored_energy_j: f64,
    pub end_stored_energy_j: f64,
    /// Direct original increments, never subtraction of rounded endpoint totals.
    pub stored_energy_change_j: f64,
    pub kelvin_heat_j: f64,
    pub maxwell_heat_j: f64,
    pub material_heat_j: f64,
    pub stress_work_j: f64,
    pub energy_balance_residual_j: f64,
}

/// Complete unpublished candidate. A refused or cancelled interval changes no
/// accepted displacement, source clock, branch history or cumulative account.
pub struct TensorShellMaterialCausalStep<'a> {
    pub history: TensorShellMaterialCausalHistory<'a>,
    pub field: TensorShellMaterialCausalField,
}

fn finite_clock(time: Dd) -> bool {
    time.hi.is_finite() && time.lo.is_finite()
}

fn material_error(
    source: &TensorPronyShell,
    element: usize,
    layer: usize,
    error: TensorPronyError,
) -> TensorShellHarmonicError {
    match error {
        TensorPronyError::Cancelled => PlateError::Cancelled.into(),
        TensorPronyError::Invalid(what) => {
            source.sections[element].material_error(element, layer, ViscoError::Parameters { what })
        }
    }
}

impl TensorPronyShell {
    /// Prepare source-bound causal plane-stress evolution of every original layer.
    /// Original CST/DKT rows, layer axes, signed thickness/reference surfaces and
    /// stationary support maps are retained. This controller admits no nonzero
    /// Kelvin viscosity. Material parameters and temperature remain fixed.
    ///
    /// Payload covers this observer, one accepted and one candidate twelve-point
    /// bank per layer, the returned six-point field, and peak temporal/query
    /// scratch. Caller input phasors and the original shell are excluded. Dense
    /// preparation storage is charged at each control cap; temporal storage is
    /// charged once at its peak cap because layers are evaluated sequentially.
    /// Work conservatively charges control.max_preparation_work plus thirteen
    /// temporal.maximum_scalar_work per layer (one preparation and twelve point
    /// queries), plus 32768*(F+1) per layer, 65536*(F+1) per facet and
    /// 256*(F+1) per full coordinate for original field recovery and receipts.
    /// Declared sub-caps, rather than guessed material ranks, bound all work.
    pub fn prepare_material_causal_observer(
        &self,
        budget: TensorShellMaterialCausalBudget,
        cx: &Cx<'_>,
    ) -> Result<TensorShellMaterialCausalObserver<'_>, TensorShellHarmonicError> {
        checkpoint(cx)?;
        if budget.maximum_frequencies > MAX_HARMONIC_MATERIAL_LINES {
            return Err(invalid(
                "causal shell material supports at most 64 authored carriers",
            ));
        }
        let f = budget.maximum_frequencies as u128;
        let facets = self.sections.len() as u128;
        let mut points = 0_u128;
        let mut bytes = core::mem::size_of::<TensorShellMaterialCausalObserver<'_>>() as u128
            + 2 * core::mem::size_of::<TensorShellMaterialCausalHistory<'_>>() as u128
            + core::mem::size_of::<TensorShellMaterialCausalField>() as u128
            + 16 * self.model.dof_map.len() as u128
            + facets
                * (core::mem::size_of::<PreparedFacet>()
                    + core::mem::size_of::<Vec<LayerHistory>>()
                    + 2 * core::mem::size_of::<Vec<[TensorPronyHistory; 12]>>())
                    as u128
            + f * (core::mem::size_of::<HarmonicPhaseInterval>()
                + 2 * core::mem::size_of::<f64>()
                + core::mem::size_of::<usize>()
                + core::mem::size_of::<[[C64; 6]; 3]>()
                + 2 * core::mem::size_of::<[C64; 3]>()) as u128
            + budget.temporal.maximum_payload_bytes as u128;
        let mut work = (f + 1) * (65536 * facets + 256 * self.model.dof_map.len() as u128);
        for section in &self.sections {
            checkpoint(cx)?;
            for layer in 0..section.layer_count() {
                checkpoint(cx)?;
                let branches = section.layer_law(layer)?.branches().len() as u128;
                points += 6;
                // Original model clone/factors, 1 virgin + 12 old + 12 new + 1
                // temporary query history. No numerical branch is discarded.
                bytes += core::mem::size_of::<LayerHistory>() as u128
                    + 4096 * (branches + 1)
                    + 8 * budget.control.max_dense_scalars as u128
                    + 26 * (core::mem::size_of::<TensorPronyHistory>() as u128 + 48 * branches)
                    + 6 * (core::mem::size_of::<TensorShellMaterialCausalPoint>() as u128
                        + f * core::mem::size_of::<[C64; 3]>() as u128);
                work += budget.control.max_preparation_work as u128
                    + 13 * budget.temporal.maximum_scalar_work as u128
                    + 32768 * (f + 1);
            }
        }
        cap(
            "causal shell material points",
            points,
            budget.field.maximum_points as u128,
        )?;
        cap(
            "causal shell material payload bytes",
            bytes,
            budget.field.maximum_payload_bytes as u128,
        )?;
        cap(
            "causal shell material work",
            work,
            budget.field.maximum_work as u128,
        )?;
        let mut prepared_facets = Vec::with_capacity(self.sections.len());
        let mut layers = Vec::with_capacity(self.sections.len());
        for (element, section) in self.sections.iter().enumerate() {
            checkpoint(cx)?;
            prepared_facets.push(PreparedFacet::new(self, element)?);
            let mut prepared_layers = Vec::with_capacity(section.layer_count());
            for layer in 0..section.layer_count() {
                checkpoint(cx)?;
                let (faces_m, to_principal) = section.layer_geometry(element, layer)?;
                let virgin = section.layer_law(layer)?.history();
                let control = PreparedTensorMixedControl::new(
                    &virgin,
                    [false, false, true, false, true, true],
                    budget.control,
                    || cx.checkpoint().is_err(),
                )
                .map_err(|e| material_error(self, element, layer, e))?;
                prepared_layers.push(LayerHistory {
                    control,
                    virgin,
                    faces_m,
                    to_principal,
                });
            }
            layers.push(prepared_layers);
        }
        checkpoint(cx)?;
        Ok(TensorShellMaterialCausalObserver {
            source: self,
            binding: Arc::new(()),
            facets: prepared_facets,
            layers,
            points: points as usize,
            budget,
        })
    }
}

impl<'a> TensorShellMaterialCausalObserver<'a> {
    pub const fn source(&self) -> &'a TensorPronyShell {
        self.source
    }

    /// Start the undeformed, relaxed original shell at a specified source clock.
    /// There is no inferred preload or instantaneous jump. The first prescribed
    /// nodal field must join zero, just as later intervals must join their history.
    pub fn start_relaxed(
        &self,
        start: Dd,
        cx: &Cx<'_>,
    ) -> Result<TensorShellMaterialCausalHistory<'a>, TensorShellHarmonicError> {
        checkpoint(cx)?;
        let start = Dd::from_pair(start.hi, start.lo);
        if !finite_clock(start) {
            return Err(invalid(
                "finite normalized causal shell source clock required",
            ));
        }
        let mut histories = Vec::with_capacity(self.layers.len());
        for layers in &self.layers {
            checkpoint(cx)?;
            let mut states = Vec::with_capacity(layers.len());
            for layer in layers {
                checkpoint(cx)?;
                states.push(core::array::from_fn(|_| layer.virgin.clone()));
            }
            histories.push(states);
        }
        checkpoint(cx)?;
        Ok(TensorShellMaterialCausalHistory {
            source: self.source,
            binding: Arc::clone(&self.binding),
            time_s: start,
            full_displacement: vec![0.0; self.source.model.dof_map.len()],
            histories,
            stored_j: 0.0,
            work_j: 0.0,
            dissipated_j: 0.0,
        })
    }

    /// Advance a continuous prescribed coherent field under exp(-i omega t).
    /// Each tone and the static bias use the source's original FREE translations
    /// and axial rotations; supports remain stationary. The complete full nodal
    /// field must join the accepted history at its retained DD instant. An empty
    /// bank holds the supplied displacement while original memory relaxes.
    ///
    /// The original coupled material controller evolves all transverse strains
    /// at zero transverse traction. Six vertex/face states bound the whole layer
    /// throughout the interval; six positive triangle/thickness quadrature states
    /// independently integrate its quadratic physical work, storage and heat.
    /// No kinetic/drilling/boundary energy or dynamic force equilibrium is inferred.
    /// Any late refusal or cancellation discards the complete candidate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_window(
        &self,
        accepted: &TensorShellMaterialCausalHistory<'a>,
        angular_frequencies_rad_s: &[f64],
        displacements: &[&[C64]],
        bias: &[f64],
        duration_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TensorShellMaterialCausalStep<'a>, TensorShellHarmonicError> {
        checkpoint(cx)?;
        if !core::ptr::eq(self.source, accepted.source)
            || !Arc::ptr_eq(&self.binding, &accepted.binding)
        {
            return Err(invalid(
                "causal shell history belongs to a different original observer",
            ));
        }
        let end_time = accepted.time_s + Dd::from_f64(duration_s);
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !finite_clock(end_time)
            || !accepted.time_s.lt(end_time)
            || angular_frequencies_rad_s.len() > self.budget.maximum_frequencies
            || angular_frequencies_rad_s
                .iter()
                .any(|w| !w.is_finite() || *w <= 0.0)
            || displacements.len() != angular_frequencies_rad_s.len()
            || bias.len() != self.source.model.free
            || bias.iter().any(|v| !v.is_finite())
            || displacements.iter().any(|u| {
                u.len() != self.source.model.free
                    || u.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
            })
        {
            return Err(invalid(
                "causal shell needs finite original free fields and a positive advancing DD interval",
            ));
        }
        let mut distinct = Vec::with_capacity(angular_frequencies_rad_s.len());
        let mut groups = Vec::with_capacity(angular_frequencies_rad_s.len());
        let mut phases = Vec::with_capacity(angular_frequencies_rad_s.len());
        for &omega in angular_frequencies_rad_s {
            let group = if let Some(i) = distinct.iter().position(|&w| w == omega) {
                i
            } else {
                let i = distinct.len();
                distinct.push(omega);
                phases.push(
                    phase_angular_interval_from_clock_dd(omega, accepted.time_s, duration_s)
                        .map_err(|_| invalid("causal shell nodal phase is unrepresentable"))?,
                );
                i
            };
            groups.push(group);
        }
        let mut endpoint = Vec::with_capacity(self.source.model.dof_map.len());
        for (full, &free) in self.source.model.dof_map.iter().enumerate() {
            checkpoint(cx)?;
            let Some(free) = free else {
                endpoint.push(0.0);
                continue;
            };
            // Coalesce carriers BEFORE endpoint evaluation and tolerance scaling:
            // a huge cancelled drive must not hide a real nodal discontinuity.
            let mut coherent = [[Dd::ZERO; 2]; MAX_HARMONIC_MATERIAL_LINES];
            for (tone, &group) in displacements.iter().zip(&groups) {
                coherent[group][0] = coherent[group][0] + Dd::from_f64(tone[free].re);
                coherent[group][1] = coherent[group][1] + Dd::from_f64(tone[free].im);
            }
            let mut start = Dd::ZERO;
            let mut end = Dd::ZERO;
            let mut scale = bias[free].abs().max(accepted.full_displacement[full].abs());
            for (z, phase) in coherent.iter().zip(&phases) {
                if !finite_clock(z[0]) || !finite_clock(z[1]) {
                    return Err(invalid(
                        "coherent causal shell displacement is unrepresentable",
                    ));
                }
                scale = scale.max(z[0].to_f64().abs()).max(z[1].to_f64().abs());
                start = start + z[0] * phase.start[0] - z[1] * phase.start[1];
                end = end + z[0] * phase.end[0] - z[1] * phase.end[1];
            }
            start = start + Dd::from_f64(bias[free]);
            end = end + Dd::from_f64(bias[free]);
            let difference = (start - Dd::from_f64(accepted.full_displacement[full]))
                .to_f64()
                .abs();
            let tolerance = 128.0 * f64::EPSILON * (2 * phases.len() + 2) as f64;
            if !finite_clock(start)
                || !finite_clock(end)
                || !difference.is_finite()
                || (scale == 0.0 && difference != 0.0)
                || (scale > 0.0 && difference / scale > tolerance)
            {
                return Err(invalid(
                    "full nodal displacement is discontinuous with its causal shell history",
                ));
            }
            endpoint.push(end.to_f64());
        }
        let mut points = Vec::with_capacity(self.points);
        let mut histories = Vec::with_capacity(self.layers.len());
        // Interval start/end storage, direct change, Kelvin/Maxwell/total heat,
        // independent stress work, and original cumulative point work/heat.
        let mut totals = [Dd::ZERO; 9];
        for (element, layers) in self.layers.iter().enumerate() {
            checkpoint(cx)?;
            let facet = &self.facets[element];
            let midpoint: Vec<_> = displacements.iter().map(|u| facet.strains(u)).collect();
            let bias_midpoint = facet.strains_by(|i| C64::from_re(bias[i]));
            let mut layer_histories = Vec::with_capacity(layers.len());
            for (layer, original) in layers.iter().enumerate() {
                checkpoint(cx)?;
                let window = original
                    .control
                    .prepare_plane_stress_harmonic_window(
                        MaterialPlane::Axes12,
                        angular_frequencies_rad_s,
                        accepted.time_s,
                        duration_s,
                        self.budget.temporal,
                        || cx.checkpoint().is_err(),
                    )
                    .map_err(|e| material_error(self.source, element, layer, e))?;
                let [bottom, top] = original.faces_m;
                let half = 0.5 * (top - bottom);
                let center = bottom + half;
                let volume = (facet.area_m2 / 3.0) * half;
                if !half.is_finite() || half <= 0.0 || !volume.is_finite() || volume <= 0.0 {
                    return Err(invalid(
                        "original causal shell quadrature volume is unrepresentable",
                    ));
                }
                let mut next = Vec::with_capacity(12);
                for index in 0..12 {
                    checkpoint(cx)?;
                    let physical = index >= 6;
                    let triangle = if physical { (index - 6) / 2 } else { index / 2 };
                    let thickness = index % 2;
                    let z = if physical {
                        center + [-1.0, 1.0][thickness] * GAUSS * half
                    } else {
                        original.faces_m[thickness]
                    };
                    let at = |rows: &[[C64; 6]; 3]| -> [C64; 3] {
                        let strain = if physical {
                            rows[triangle]
                        } else {
                            core::array::from_fn(|i| {
                                rows[triangle][i] + rows[(triangle + 2) % 3][i]
                                    - rows[(triangle + 1) % 3][i]
                            })
                        };
                        principal_strain(&strain, z, &original.to_principal)
                    };
                    let strain: Vec<_> = midpoint.iter().map(at).collect();
                    let point_bias = at(&bias_midpoint).map(|z| z.re);
                    let response = window
                        .evaluate(
                            &strain,
                            &point_bias,
                            &accepted.histories[element][layer][index],
                            || cx.checkpoint().is_err(),
                        )
                        .map_err(|e| material_error(self.source, element, layer, e))?;
                    if physical {
                        let material = response.material;
                        let values = [
                            material.start_stored_energy_j_m3,
                            material.end_stored_energy_j_m3,
                            material.stored_energy_change_j_m3,
                            material.kelvin_heat_j_m3,
                            material.maxwell_heat_j_m3,
                            material.material_heat_j_m3,
                            material.stress_work_j_m3,
                            response.history.work(),
                            response.history.dissipated(),
                        ];
                        for (sum, value) in totals.iter_mut().zip(values) {
                            *sum = *sum + Dd::from_f64(volume) * Dd::from_f64(value);
                        }
                        points.push(TensorShellMaterialCausalPoint {
                            element,
                            layer,
                            triangle_quadrature_point: triangle,
                            thickness_quadrature_point: thickness,
                            triangle_barycentric: BARYCENTRIC[triangle],
                            thickness_position_m: z,
                            volume_m3: volume,
                            in_plane_strain: strain,
                            in_plane_bias: point_bias,
                            end_strain: response.history.strain(),
                            end_stress: response.history.stress(),
                            strain_bounds: response.strain_bounds,
                            material,
                        });
                    }
                    next.push(response.history);
                }
                let bank: [TensorPronyHistory; 12] = next
                    .try_into()
                    .map_err(|_| invalid("incomplete original causal shell layer history"))?;
                layer_histories.push(bank);
            }
            histories.push(layer_histories);
        }
        let values = totals.map(Dd::to_f64);
        let balance = (totals[6] - totals[2] - totals[5]).to_f64();
        if values.iter().any(|v| !v.is_finite())
            || !balance.is_finite()
            || [0, 1, 3, 4, 5, 8].iter().any(|&i| values[i] < 0.0)
        {
            return Err(invalid("causal shell material account is unrepresentable"));
        }
        checkpoint(cx)?;
        Ok(TensorShellMaterialCausalStep {
            history: TensorShellMaterialCausalHistory {
                source: self.source,
                binding: Arc::clone(&self.binding),
                time_s: end_time,
                full_displacement: endpoint,
                histories,
                stored_j: values[1],
                work_j: values[7],
                dissipated_j: values[8],
            },
            field: TensorShellMaterialCausalField {
                start_time_s_dd: accepted.time_s,
                end_time_s_dd: end_time,
                duration_s,
                angular_frequencies_rad_s: angular_frequencies_rad_s.to_vec(),
                points,
                start_stored_energy_j: values[0],
                end_stored_energy_j: values[1],
                stored_energy_change_j: values[2],
                kelvin_heat_j: values[3],
                maxwell_heat_j: values[4],
                material_heat_j: values[5],
                stress_work_j: values[6],
                energy_balance_residual_j: balance,
            },
        })
    }
}
