//! Explicit physical inputs for the two-interval measurement-fitting consumer.
//!
//! This is a bounded example input format, not FrankenScript or a material-card
//! authority. Geometry, conductivity, caloric data and exterior heat loss are
//! caller declarations; all numerical assembly and evolution stay in fs crates.
use std::collections::{BTreeMap, BTreeSet};
use fs_blake3::ContentHash;
use fs_conduction::{ConductionMesh, ConductivityModel, ThermalBc, COMPONENT_POWER_DIMS};
use fs_conduction::bc::history::march::moving::{MovingHeatBudget, MovingHeatSource};
use fs_exec::Cx;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_qty::{Dims, QtyAny};
use fs_rep_mesh::TetComplex;
use fs_scenario::signal::{Interp, TimeSignal};

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

/// Input-shape bounds, not a total-memory or physical-accuracy certificate.
#[derive(Debug, Clone, Copy)]
pub struct SpecimenBudget {
    pub max_bytes: usize,
    pub max_vertices: usize,
    pub max_elements: usize,
    pub max_phase_knots: usize,
}

/// Owned inputs. Borrow these same objects for every candidate and its replay.
pub struct ScanSpecimen {
    pub mesh: ConductionMesh,
    pub material: ConductivityModel,
    pub curve: EquilibriumEnthalpyPhaseCurve,
    pub exterior: ThermalBc,
    pub reference_density_kg_m3: f64,
    pub initial_specific_enthalpy_j_kg: f64,
    pub switch_time_s: f64,
    pub end_time_s: f64,
    pub first_center_m: [f64; 3],
    pub second_center_m: [f64; 3],
    pub half_width_m: [f64; 3],
    /// Identity of supplied bytes, not proof of material calibration or geometry.
    pub identity: ContentHash,
}

impl ScanSpecimen {
    /// Admit every field before constructing the original mesh/material objects.
    /// Repeated scalar declarations, unknown keys and implicit defaults refuse.
    /// The topology must be a supplied conforming tetrahedral mesh. This loader
    /// rejects orphan vertices, duplicate tets, nonmanifold faces and degenerate
    /// elements; it does not detect intersecting/overlapping unrelated tets.
    pub fn from_text(cx: &Cx<'_>, text: &str, budget: SpecimenBudget) -> Result<Self> {
        poll(cx)?;
        if text.len() > budget.max_bytes { return Err("specimen exceeds its byte allowance".into()); }
        let mut lines = text.lines().enumerate().filter_map(|(i, line)| {
            let line = line.split('#').next().unwrap_or("").trim();
            (!line.is_empty()).then_some((i + 1, line))
        });
        if lines.next().map(|(_, line)| line) != Some("FS_SCAN_SPECIMEN_V1") {
            return Err("specimen must begin with FS_SCAN_SPECIMEN_V1".into());
        }
        let mut seen = BTreeSet::new();
        let (mut conductivity, mut density, mut initial, mut switch, mut end) = (None, None, None, None, None);
        let (mut first_center, mut second_center, mut half_width, mut exterior) = (None, None, None, None);
        let mut positions = Vec::new();
        let mut tets = Vec::new();
        let mut knots = Vec::<[f64; 3]>::new();
        for (line_number, line) in lines {
            poll(cx)?;
            let result = (|| -> Result<()> {
                let fields: Vec<_> = line.split_whitespace().take(6).collect();
                let key = fields[0];
                if !matches!(key, "vertex_m" | "tet" | "phase_knot") && !seen.insert(key) {
                    return Err(format!("duplicate specimen field {key}").into());
                }
                let args = &fields[1..];
                match key {
                    "conductivity_w_m_k" => conductivity = Some(numbers::<1>(args)?[0]),
                    "reference_density_kg_m3" => density = Some(numbers::<1>(args)?[0]),
                    "initial_specific_enthalpy_j_kg" => initial = Some(numbers::<1>(args)?[0]),
                    "switch_time_s" => switch = Some(numbers::<1>(args)?[0]),
                    "end_time_s" => end = Some(numbers::<1>(args)?[0]),
                    "first_center_m" => first_center = Some(numbers::<3>(args)?),
                    "second_center_m" => second_center = Some(numbers::<3>(args)?),
                    "half_width_m" => half_width = Some(numbers::<3>(args)?),
                    "exterior" => exterior = Some(match args {
                        ["adiabatic"] => ThermalBc::adiabatic(),
                        ["robin", htc, ambient] => {
                            let [h, t] = numbers::<2>(&[htc, ambient])?;
                            if h <= 0.0 || t <= 0.0 { return Err("Robin h and ambient kelvin must be positive".into()); }
                            ThermalBc::robin(h, t)?
                        }
                        _ => return Err("exterior is adiabatic or robin HTC_W_m2_K AMBIENT_K".into()),
                    }),
                    "vertex_m" => {
                        if positions.len() >= budget.max_vertices { return Err("specimen vertex allowance exceeded".into()); }
                        positions.push(numbers::<3>(args)?);
                    }
                    "tet" => {
                        if tets.len() >= budget.max_elements || args.len() != 4 {
                            return Err("tet needs four vertex indices within the element allowance".into());
                        }
                        tets.push([args[0].parse::<u32>()?, args[1].parse()?, args[2].parse()?, args[3].parse()?]);
                    }
                    "phase_knot" => {
                        if knots.len() >= budget.max_phase_knots { return Err("specimen phase-knot allowance exceeded".into()); }
                        knots.push(numbers::<3>(args)?);
                    }
                    _ => return Err(format!("unknown specimen field {key}").into()),
                }
                Ok(())
            })();
            result.map_err(|error| format!("specimen line {line_number}: {error}"))?;
        }
        let conductivity = conductivity.ok_or("missing conductivity_w_m_k")?;
        let density = density.ok_or("missing reference_density_kg_m3")?;
        let initial = initial.ok_or("missing initial_specific_enthalpy_j_kg")?;
        let switch = switch.ok_or("missing switch_time_s")?;
        let end = end.ok_or("missing end_time_s")?;
        let first_center = first_center.ok_or("missing first_center_m")?;
        let second_center = second_center.ok_or("missing second_center_m")?;
        let half_width: [f64; 3] = half_width.ok_or("missing half_width_m")?;
        let exterior = exterior.ok_or("missing exterior")?;
        if conductivity <= 0.0 || density <= 0.0 || switch <= 0.0 || end <= switch
            || !(end - switch).is_finite() || half_width.iter().any(|&x| x <= 0.0)
        { return Err("positive conductivity/density/widths and 0 < switch < end are required".into()); }
        if positions.len() < 4 || tets.is_empty() || knots.len() < 2 {
            return Err("specimen needs four vertices, a tetrahedron and two caloric knots".into());
        }
        if initial < knots[0][0] || initial > knots[knots.len()-1][0] {
            return Err("initial enthalpy is outside the declared caloric interval".into());
        }
        let mut used = vec![false; positions.len()];
        let mut unique = BTreeSet::new();
        let mut faces = BTreeMap::<[u32; 3], usize>::new();
        for tet in &tets {
            poll(cx)?;
            let mut key = *tet;
            key.sort_unstable();
            if key.windows(2).any(|pair| pair[0] == pair[1]) || !unique.insert(key) {
                return Err("tetrahedra must have four distinct vertices and unique connectivity".into());
            }
            for omitted in 0..4 {
                let mut face = [0; 3];
                let mut at = 0;
                for (i, &vertex) in key.iter().enumerate() {
                    if i != omitted { face[at] = vertex; at += 1; }
                }
                let count = faces.entry(face).or_default();
                *count += 1;
                if *count > 2 { return Err("a mesh face belongs to more than two tetrahedra".into()); }
            }
            for &vertex in tet {
                *used.get_mut(vertex as usize).ok_or("tet vertex is outside supplied positions")? = true;
            }
        }
        if used.contains(&false) { return Err("specimen has an unused vertex with no thermal mass".into()); }
        let identity = fs_blake3::hash_domain("scan-fit-declared-specimen-v1", text.as_bytes());
        let phase_knots = knots.into_iter().map(|[h, t, f]| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h, temperature_k: t, liquid_mass_fraction: f,
            bulk_density_kg_m3: density,
        }).collect();
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(identity, phase_knots)?;
        let material = ConductivityModel::isotropic_declared(conductivity)?;
        poll(cx)?;
        // Original constructors own element geometry, boundary extraction and
        // constitutive checks. Construction work is bounded by the input caps.
        let mesh = ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions)?;
        poll(cx)?;
        Ok(Self { mesh, material, curve, exterior, reference_density_kg_m3: density,
            initial_specific_enthalpy_j_kg: initial, switch_time_s: switch, end_time_s: end,
            first_center_m: first_center, second_center_m: second_center,
            half_width_m: half_width, identity })
    }

    /// The historical manufactured specimen, retained for existing commands.
    pub fn synthetic_reference(cx: &Cx<'_>) -> Result<Self> {
        poll(cx)?;
        let (complex, positions) = fs_conduction::fixtures::box_grid([1, 1, 1], [1.0; 3]);
        let density = 1.0;
        let identity = ContentHash([76; 32]);
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(identity,
            [(-200.0, 280.0, 0.0), (100.0, 310.0, 0.0), (200.0, 310.0, 1.0), (900.0, 380.0, 1.0)]
                .into_iter().map(|(h, t, f)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: h, temperature_k: t,
                    liquid_mass_fraction: f, bulk_density_kg_m3: density,
                }).collect())?;
        Ok(Self { mesh: ConductionMesh::new(complex, positions)?,
            material: ConductivityModel::isotropic_declared(2.0)?, curve,
            exterior: ThermalBc::adiabatic(), reference_density_kg_m3: density,
            initial_specific_enthalpy_j_kg: 150.0, switch_time_s: 1.0, end_time_s: 2.0,
            first_center_m: [0.25, 0.5, 0.5], second_center_m: [0.75, 0.5, 0.5],
            half_width_m: [1.4, 1.2, 1.1], identity })
    }

    /// Original held schedule: two fitted absorbed powers, then exact zero.
    /// Position also switches explicitly; no finite-speed claim is made.
    pub fn source(&self, initial_w: [f64; 2]) -> Result<MovingHeatSource> {
        if initial_w.iter().any(|x| !x.is_finite() || *x < 0.0) {
            return Err("initial absorbed powers must be finite nonnegative watts".into());
        }
        let times = vec![0.0, self.switch_time_s, self.end_time_s];
        let length = Dims([1, 0, 0, 0, 0, 0]);
        let center = std::array::from_fn(|axis| {
            if self.first_center_m[axis].to_bits() == self.second_center_m[axis].to_bits() {
                TimeSignal::Constant(QtyAny::new(self.first_center_m[axis], length))
            } else {
                TimeSignal::Table { times: times.clone(),
                    values: vec![self.first_center_m[axis], self.second_center_m[axis], self.second_center_m[axis]],
                    dims: length, interp: Interp::Hold }
            }
        });
        Ok(MovingHeatSource::new(center, self.half_width_m,
            TimeSignal::Table { times, values: vec![initial_w[0], initial_w[1], 0.0],
                dims: COMPONENT_POWER_DIMS, interp: Interp::Hold }, 0.0,
            MovingHeatBudget { max_vertices: self.mesh.vertex_count(), max_elements: self.mesh.element_count(),
                min_active_vertices: 1, power_relative_tolerance: 1e-11 })?)
    }
}

fn numbers<const N: usize>(args: &[&str]) -> Result<[f64; N]> {
    if args.len() != N { return Err(format!("expected {N} numeric values").into()); }
    let mut out = [0.0; N];
    for (dst, text) in out.iter_mut().zip(args) {
        *dst = text.parse()?;
        if !dst.is_finite() { return Err("specimen numbers must be finite SI values".into()); }
    }
    Ok(out)
}
fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| "specimen loading cancelled".into())
}

#[cfg(test)]
mod tests;
