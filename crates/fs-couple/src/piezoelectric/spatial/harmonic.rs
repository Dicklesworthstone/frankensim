//! Spatial piezoelectric frequency response with physical electrode shunts.
//!
//! Peak phasors use `exp(-i omega t)`, as in the native structural harmonic
//! owners. The SAME prepared K, M, G and K_e solve
//! `[K-omega^2 M, G; G^T, -K_e-i Y/omega] [u,phi] = [f,-i I/omega]`.
//! Y is assembled from source-bound external R, L and C elements. There is no
//! modal truncation, fitted actuator gain or invented ceramic damping.
//!
//! A complex symmetric pencil is represented by the real symmetric block
//! `[Re(A),-Im(A);-Im(A),-Re(A)]` and uses the existing sparse MINRES owner.
//! Scaling uses positive stiffness/inertia magnitudes, so cancellation of a
//! dynamic diagonal is not mistaken for a singular system. Physical field and
//! charge residuals and mean work are independently reconstructed afterwards.
//!
//! This is a periodic particular solution of the retained, fixed-material mesh.
//! It does not establish startup decay, uniqueness at an unexcited lossless
//! pole, mesh convergence, intrinsic coupling loss, radiation or empirical
//! material accuracy. A residual is not a forward-error bound near resonance.
//!
//! [`TetPiezoelectricModel::solve_harmonic_with_relaxation`] additionally admits
//! original-card Maxwell mechanical and directed Debye polarization processes.
//! Their complex increments use the same mesh gradients and reciprocal G;
//! independently recovered material memory and heat extend the original work
//! account. The fixed real piezoelectric coupling itself remains reversible.

use super::transient::load::{ResistiveElectrodeLoad, reactive::ReactiveElectrodeLoad};
use super::{
    BoundaryMap, ElectrodeCondition, PiezoelectricElectrode, Result, TetPiezoelectricControls,
    TetPiezoelectricModel, equilibrated_solve_scaled, invalid, poll, positive, prepare_boundary,
    validate_controls,
};
use fs_exec::Cx;
use fs_material::piezoelectric::tensor::relaxation::RelaxingTensorPiezoelectric;
use fs_math::{c64::C64, det};
use fs_sparse::Coo;

mod relaxation;

/// Prescribed peak phasor. A zero current is an AC open circuit, and a zero
/// voltage is a short circuit to the component's explicit voltage reference.
#[derive(Clone, Copy, Debug)]
pub enum HarmonicElectrodeCondition {
    /// Electrode potential [V].
    Voltage(C64),
    /// External source current INTO the electrode [A], before shunt outflow.
    Current(C64),
}

/// Equipotential complete exterior faces, under the original spatial rules.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicElectrode<'a> {
    /// Original mesh vertices forming complete, disjoint exterior faces.
    pub nodes: &'a [usize],
    /// Prescribed potential or external source current.
    pub condition: HarmonicElectrodeCondition,
}

/// Original source-bound circuit element, evaluated at the requested frequency.
#[derive(Clone, Debug)]
pub enum HarmonicElectrodeLoad {
    /// Uniform fixed-condition Ohmic conductor.
    Resistive(ResistiveElectrodeLoad),
    /// Original material winding or external dielectric capacitor.
    Reactive(ReactiveElectrodeLoad),
}

impl HarmonicElectrodeLoad {
    /// Current is positive from the first electrode to the second.
    pub fn electrodes(&self) -> [usize; 2] {
        match self {
            Self::Resistive(load) => load.electrodes,
            Self::Reactive(load) => load.electrodes(),
        }
    }

    /// Source-derived admittance [S] under `exp(-i omega t)`. Source conditions
    /// stay fixed; an arbitrary complex gain cannot masquerade as a material.
    pub fn admittance_s(&self, omega_rad_s: f64) -> Result<C64> {
        if !positive(omega_rad_s) {
            return Err(invalid(
                "harmonic loads require positive finite angular frequency",
            ));
        }
        let resistance = |conductor: &fs_material::conductor::ResolvedConductor| {
            if !conductor
                .source_property()
                .is_static_at_fixed_external_conditions()
            {
                return Err(invalid("harmonic shunt requires a fixed Ohmic source"));
            }
            Ok(conductor.resistance_ohm())
        };
        let value = match self {
            Self::Resistive(load) => C64::from_re(1.0 / resistance(&load.conductor)?),
            Self::Reactive(ReactiveElectrodeLoad::SeriesRl {
                inductor,
                conductor,
                ..
            }) => {
                let r = match conductor {
                    Some(source) => resistance(source)?,
                    None => 0.0,
                };
                let reactance = omega_rad_s * inductor.inductance_h();
                if !positive(reactance) {
                    return Err(invalid("inductive reactance is not representable"));
                }
                C64::new(r, -reactance).recip()
            }
            Self::Reactive(ReactiveElectrodeLoad::Capacitor { capacitor, .. }) => {
                let susceptance = omega_rad_s * capacitor.capacitance_f();
                if !positive(susceptance) {
                    return Err(invalid("capacitive susceptance is not representable"));
                }
                C64::new(0.0, -susceptance)
            }
        };
        if !finite(value) || value.abs() == 0.0 || value.re < 0.0 {
            return Err(invalid("passive source admittance is not representable"));
        }
        Ok(value)
    }
}

/// Sparse solve and full-cycle physical admission. No dense matrix is allocated.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicPiezoelectricControls {
    /// Original residual, iteration, storage and small-signal field controls.
    pub field: TetPiezoelectricControls,
    /// Absolute cycle-averaged supply minus load heat tolerance [W].
    pub power_tolerance_w: f64,
    /// Maximum doubled real coordinate count, checked before matrix assembly.
    pub maximum_real_coordinates: usize,
}

/// Complete peak phasors in the original mesh/material frame.
#[derive(Clone, Debug)]
pub struct HarmonicPiezoelectricCell {
    /// Peak infinitesimal strain in the original Mandel frame.
    pub strain_mandel: [C64; 6],
    /// Peak work-conjugate tensile stress [Pa].
    pub stress_mandel_pa: [C64; 6],
    /// Peak electric field, equal to minus the potential gradient [V/m].
    pub electric_field_v_m: [C64; 3],
    /// Peak electric displacement from the reciprocal source tensor [C/m²].
    pub displacement_c_m2: [C64; 3],
    /// Mean physical Maxwell heat integrated over this original cell [W].
    pub average_mechanical_heat_w: f64,
    /// Mean physical Debye polarization heat over this original cell [W].
    pub average_dielectric_heat_w: f64,
    /// Mean mechanical plus polarization memory storage in this cell [J].
    /// Already included in the response's total body stored energy.
    pub average_memory_energy_j: f64,
}

/// Body free charge and external source current, including attached load current.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicElectrodeResponse {
    /// Absolute prescribed or solved peak potential [V].
    pub voltage_v: C64,
    /// Total body free charge, excluding separate volume charge [C].
    pub charge_c: C64,
    /// External current into the electrode, before shunt outflow [A].
    pub source_current_a: C64,
}

/// Source, actual load flow and mean heat. Reactive storage does not produce heat.
#[derive(Clone, Debug)]
pub struct HarmonicElectrodeLoadResponse {
    /// Original immutable physical source and electrode orientation.
    pub load: HarmonicElectrodeLoad,
    /// Peak voltage from first to second electrode [V].
    pub voltage_v: C64,
    /// Peak current from first to second electrode [A].
    pub current_a: C64,
    /// Mean Joule heat in this branch's original conductor [W].
    pub average_heat_w: f64,
}

/// Periodic displacement, electric fields and independently checked mean work.
#[derive(Clone, Debug)]
pub struct HarmonicPiezoelectricResponse {
    /// Peak original-node displacements [m].
    pub displacement_m: Vec<[C64; 3]>,
    /// Peak velocities, equal to -i omega times displacement [m/s].
    pub velocity_m_s: Vec<[C64; 3]>,
    /// Absolute peak node potentials [V]; work uses retained differences.
    pub voltage_v: Vec<C64>,
    /// Dynamic internal force minus applied force [N]; free entries are residuals.
    pub support_force_n: Vec<[C64; 3]>,
    /// Peak constitutive state in original tetrahedron order.
    pub cells: Vec<HarmonicPiezoelectricCell>,
    /// Potentials, body charge and source current in original electrode order.
    pub electrodes: Vec<HarmonicElectrodeResponse>,
    /// Actual branch current and heat in original load order.
    pub loads: Vec<HarmonicElectrodeLoadResponse>,
    /// Cycle-averaged elastic, electric and material-memory body storage [J].
    pub average_stored_energy_j: f64,
    /// Cycle-averaged kinetic storage from the original consistent mass [J].
    pub average_kinetic_energy_j: f64,
    /// Applied forces AND reactions at prescribed moving supports [W].
    pub average_mechanical_power_w: f64,
    /// Electrode sources AND prescribed nodal-charge sources [W].
    pub average_electrical_power_w: f64,
    /// External resistor heat, counted once, never attributed to the ceramic [W].
    pub average_load_heat_w: f64,
    /// Intrinsic mechanical plus dielectric heat, integrated over all cells [W].
    pub average_material_heat_w: f64,
    /// Mechanical plus electrical supply minus load AND material heat [W].
    pub power_balance_residual_w: f64,
    /// Largest peak complex force residual on an unconstrained mechanical DOF [N].
    pub maximum_force_residual_n: f64,
    /// Largest peak free-node or total current-electrode charge residual [C].
    pub maximum_charge_residual_c: f64,
    /// Total native MINRES iterations, including bounded residual corrections.
    pub iterations: usize,
    /// Independently recomputed residual of the equilibrated real block system.
    pub relative_residual: f64,
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn current(q: C64, omega: f64) -> C64 {
    C64::new(omega * q.im, -omega * q.re)
}
fn charge(i: C64, omega: f64) -> C64 {
    C64::new(-i.im / omega, i.re / omega)
}
fn power(e: C64, f: C64) -> f64 {
    0.5 * (e * f.conj()).re
}
fn combine<const N: usize>(a: [f64; N], b: [f64; N]) -> [C64; N] {
    core::array::from_fn(|i| C64::new(a[i], b[i]))
}

// Largest Euclidean norm over the ENTIRE sinusoidal cycle. Checking real and
// imaginary components separately misses a sqrt(2) overrun for equal phases.
fn cycle_norm<const N: usize>(v: &[C64; N]) -> f64 {
    let scale = v
        .iter()
        .map(|z| z.re.abs().max(z.im.abs()))
        .fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return 0.0;
    }
    let (mut a, mut b, mut c) = (0.0, 0.0, 0.0);
    for z in v {
        let (x, y) = (z.re / scale, z.im / scale);
        a += x * x;
        b += y * y;
        c += x * y;
    }
    scale * det::sqrt(0.5 * (a + b + det::hypot(a - b, 2.0 * c)))
}

impl TetPiezoelectricModel<'_> {
    /// Solve all mechanical/electrical phasors at one positive angular frequency.
    /// Prescribed displacements may move harmonically; all others are free.
    /// Nodal charge is an integrated volume-charge phasor [C]. Every connected
    /// component needs a voltage-reference electrode. Loads cannot bridge
    /// distinct reference components. Both quadratures and their full-cycle
    /// strain/field envelope are checked before any response is published.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_harmonic(
        &self,
        omega_rad_s: f64,
        prescribed_displacements: &[(usize, C64)],
        electrodes: &[HarmonicElectrode<'_>],
        loads: &[HarmonicElectrodeLoad],
        nodal_forces_n: &[[C64; 3]],
        nodal_free_charge_c: &[C64],
        controls: HarmonicPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<HarmonicPiezoelectricResponse> {
        self.solve_harmonic_internal(
            omega_rad_s,
            None,
            prescribed_displacements,
            electrodes,
            loads,
            nodal_forces_n,
            nodal_free_charge_c,
            controls,
            cx,
        )
    }

    /// Solve with source-bound intrinsic mechanical and dielectric relaxation.
    /// Supply one uniform law or one per original tetrahedron. Each law must
    /// retain the EXACT original base tensor identity and mesh frame. Its base
    /// C^E is explicitly relaxed and epsilon^S is explicitly high-frequency;
    /// no single-frequency datasheet values are reinterpreted automatically.
    ///
    /// The same sparse pencil, reference inertia, coupling, boundary reduction,
    /// work audit and physical residual gates are retained. Source frequency
    /// domains are intersected by the material owner. Heat is reconstructed from
    /// branch states, with external load heat reported separately. This is a
    /// held-temperature linear model, not thermal feedback or hysteresis.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_harmonic_with_relaxation(
        &self,
        omega_rad_s: f64,
        intrinsic_materials: &[RelaxingTensorPiezoelectric],
        prescribed_displacements: &[(usize, C64)],
        electrodes: &[HarmonicElectrode<'_>],
        loads: &[HarmonicElectrodeLoad],
        nodal_forces_n: &[[C64; 3]],
        nodal_free_charge_c: &[C64],
        controls: HarmonicPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<HarmonicPiezoelectricResponse> {
        self.solve_harmonic_internal(
            omega_rad_s,
            Some(intrinsic_materials),
            prescribed_displacements,
            electrodes,
            loads,
            nodal_forces_n,
            nodal_free_charge_c,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_harmonic_internal(
        &self,
        omega_rad_s: f64,
        intrinsic_materials: Option<&[RelaxingTensorPiezoelectric]>,
        prescribed_displacements: &[(usize, C64)],
        electrodes: &[HarmonicElectrode<'_>],
        loads: &[HarmonicElectrodeLoad],
        nodal_forces_n: &[[C64; 3]],
        nodal_free_charge_c: &[C64],
        controls: HarmonicPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<HarmonicPiezoelectricResponse> {
        poll(cx)?;
        validate_controls(controls.field)?;
        let n = self.nodes_m.len();
        let nu = 3 * n;
        let omega2 = omega_rad_s * omega_rad_s;
        if !positive(omega_rad_s)
            || !positive(omega2)
            || !positive(controls.power_tolerance_w)
            || nodal_forces_n.len() != n
            || nodal_free_charge_c.len() != n
            || nodal_forces_n
                .iter()
                .flatten()
                .chain(nodal_free_charge_c)
                .any(|&z| !finite(z))
            || prescribed_displacements.iter().any(|&(_, z)| !finite(z))
        {
            return Err(invalid(
                "finite complete phasors and positive representable frequency/power controls required",
            ));
        }
        let source_conditions = electrodes
            .iter()
            .map(|e| {
                let (voltage, value) = match e.condition {
                    HarmonicElectrodeCondition::Voltage(v) => (true, v),
                    HarmonicElectrodeCondition::Current(i) => (false, charge(i, omega_rad_s)),
                };
                if !finite(value) {
                    return Err(invalid("harmonic electrode source is not representable"));
                }
                Ok((voltage, value))
            })
            .collect::<Result<Vec<_>>>()?;
        let real = |z: C64| z.re;
        let imaginary = |z: C64| z.im;
        let component = |part: fn(C64) -> f64| -> Result<BoundaryMap> {
            let supports = prescribed_displacements
                .iter()
                .map(|&(i, z)| (i, part(z)))
                .collect::<Vec<_>>();
            let boundary = electrodes
                .iter()
                .zip(&source_conditions)
                .map(|(e, &(v, z))| PiezoelectricElectrode {
                    nodes: e.nodes,
                    condition: if v {
                        ElectrodeCondition::Voltage(part(z))
                    } else {
                        ElectrodeCondition::Charge(part(z))
                    },
                })
                .collect::<Vec<_>>();
            let forces = nodal_forces_n
                .iter()
                .map(|f| f.map(part))
                .collect::<Vec<_>>();
            let charges = nodal_free_charge_c
                .iter()
                .map(|&q| part(q))
                .collect::<Vec<_>>();
            prepare_boundary(
                self.nodes_m,
                self.tetrahedra,
                &supports,
                &boundary,
                &forces,
                &charges,
                cx,
            )
        };
        let boundary = [component(real)?, component(imaginary)?];
        let d = boundary[0].rhs.len();
        let dimensions = d
            .checked_mul(2)
            .ok_or_else(|| invalid("harmonic coordinate count overflows"))?;
        if dimensions > controls.maximum_real_coordinates {
            return Err(invalid(
                "harmonic solve exceeds its explicit real-coordinate budget",
            ));
        }
        let intrinsic_coefficients = intrinsic_materials
            .map(|materials| relaxation::prepare(self, materials, omega_rad_s, cx))
            .transpose()?;
        let mut full = boundary[0]
            .prescribed
            .iter()
            .zip(&boundary[1].prescribed)
            .map(|(&re, &im)| C64::new(re, im))
            .collect::<Vec<_>>();
        let mut rhs = boundary[0].rhs.clone();
        rhs.extend(boundary[1].rhs.iter().map(|x| -x));
        let mut matrix = Coo::new(dimensions, dimensions);
        let mut scale = vec![0.0; d];
        let mut add =
            |row: usize, column: usize, value: C64, diagonal_magnitude: f64| -> Result<()> {
                if !finite(value) || !diagonal_magnitude.is_finite() {
                    return Err(invalid("harmonic matrix coefficient is not representable"));
                }
                if let Some(r) = boundary[0].index[row] {
                    if row == column {
                        scale[r] += diagonal_magnitude;
                    }
                    if let Some(c) = boundary[0].index[column] {
                        if value.re != 0.0 {
                            matrix.push(r, c, value.re);
                            matrix.push(r + d, c + d, -value.re);
                        }
                        if value.im != 0.0 {
                            matrix.push(r, c + d, -value.im);
                            matrix.push(r + d, c, -value.im);
                        }
                    } else {
                        let shift = value * full[column];
                        rhs[r] -= shift.re;
                        rhs[r + d] += shift.im;
                    }
                }
                Ok(())
            };
        for row in 0..nu {
            poll(cx)?;
            for (operator, coefficient) in
                [(&self.solid.stiffness, 1.0), (&self.solid.mass, -omega2)]
            {
                let (cols, vals) = operator.row(row);
                for (&col, &value) in cols.iter().zip(vals) {
                    let value = value * coefficient;
                    add(row, col, C64::from_re(value), value.abs())?;
                }
            }
            let (cols, vals) = self.coupling.row(row);
            for (&col, &value) in cols.iter().zip(vals) {
                add(row, nu + col, C64::from_re(value), 0.0)?;
                add(nu + col, row, C64::from_re(value), 0.0)?;
            }
        }
        for row in 0..n {
            poll(cx)?;
            let (cols, vals) = self.dielectric.stiffness().row(row);
            for (&col, &value) in cols.iter().zip(vals) {
                add(nu + row, nu + col, C64::from_re(-value), value.abs())?;
            }
        }
        if let Some(coefficients) = &intrinsic_coefficients {
            relaxation::stamp(self, coefficients, &mut add, cx)?;
        }
        let mut admittances = Vec::with_capacity(loads.len());
        for load in loads {
            poll(cx)?;
            let [a, b] = load.electrodes();
            if a == b || a >= electrodes.len() || b >= electrodes.len() {
                return Err(invalid(
                    "harmonic loads require distinct in-range electrodes",
                ));
            }
            let (a, b) = (electrodes[a].nodes[0], electrodes[b].nodes[0]);
            if boundary[0].voltage_reference_electrode[a]
                != boundary[0].voltage_reference_electrode[b]
            {
                return Err(invalid(
                    "harmonic shunts must share one dielectric voltage-reference component",
                ));
            }
            let y = load.admittance_s(omega_rad_s)?;
            let stamp = C64::new(y.im / omega_rad_s, -y.re / omega_rad_s);
            if !finite(stamp) || stamp.abs() == 0.0 {
                return Err(invalid(
                    "harmonic load charge coefficient is not representable",
                ));
            }
            for (row, col, sign) in [(a, a, 1.0), (a, b, -1.0), (b, a, -1.0), (b, b, 1.0)] {
                add(nu + row, nu + col, stamp.scale(sign), stamp.abs())?;
            }
            admittances.push(y);
        }
        for s in &mut scale {
            *s = 1.0 / det::sqrt(*s);
        }
        scale.extend_from_within(..);
        let (solution, iterations, relative_residual) =
            equilibrated_solve_scaled(matrix.assemble(), rhs, scale, controls.field, cx)?;
        for (value, index) in full.iter_mut().zip(&boundary[0].index) {
            if let Some(i) = index {
                *value = C64::new(solution[*i], solution[*i + d]);
            }
        }
        self.recover_harmonic(
            omega_rad_s,
            &full,
            &boundary,
            electrodes,
            loads,
            &admittances,
            &source_conditions,
            nodal_forces_n,
            nodal_free_charge_c,
            intrinsic_materials,
            controls,
            iterations,
            relative_residual,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn recover_harmonic(
        &self,
        omega: f64,
        full: &[C64],
        boundary: &[BoundaryMap; 2],
        electrodes: &[HarmonicElectrode<'_>],
        loads: &[HarmonicElectrodeLoad],
        admittances: &[C64],
        source_conditions: &[(bool, C64)],
        forces: &[[C64; 3]],
        charges: &[C64],
        intrinsic_materials: Option<&[RelaxingTensorPiezoelectric]>,
        controls: HarmonicPiezoelectricControls,
        iterations: usize,
        relative_residual: f64,
        cx: &Cx<'_>,
    ) -> Result<HarmonicPiezoelectricResponse> {
        let n = self.nodes_m.len();
        let nu = 3 * n;
        let mut effective_force = forces.to_vec();
        let mut kinetic = 0.0;
        for row in 0..nu {
            poll(cx)?;
            let (cols, vals) = self.solid.mass.row(row);
            let mass_u = cols
                .iter()
                .zip(vals)
                .fold(C64::ZERO, |s, (&j, &a)| s + full[j].scale(a));
            let inertia = mass_u.scale(omega * omega);
            effective_force[row / 3][row % 3] = effective_force[row / 3][row % 3] + inertia;
            kinetic += 0.25 * (full[row].conj() * inertia).re;
        }
        let intrinsic = intrinsic_materials
            .map(|materials| relaxation::recover(self, materials, omega, full, cx))
            .transpose()?;
        // Recover the original reversible operator with the exact additional
        // constitutive force/charge moved to its right-hand side. Its residuals
        // and reactions then remain those of the COMPLETE physical equations.
        // Physical source work below still uses the original external charges.
        let mut reversible_charges = charges.to_vec();
        if let Some(intrinsic) = &intrinsic {
            for node in 0..n {
                poll(cx)?;
                for axis in 0..3 {
                    effective_force[node][axis] =
                        effective_force[node][axis] - intrinsic.force_n[node][axis];
                }
                reversible_charges[node] = reversible_charges[node] - intrinsic.charge_c[node];
            }
        }
        let mut outgoing = vec![C64::ZERO; electrodes.len()];
        let mut load_responses = Vec::with_capacity(loads.len());
        let mut heat = 0.0;
        for (load, &y) in loads.iter().zip(admittances) {
            poll(cx)?;
            let [a, b] = load.electrodes();
            let voltage = full[nu + electrodes[a].nodes[0]] - full[nu + electrodes[b].nodes[0]];
            let flow = y * voltage;
            outgoing[a] = outgoing[a] + flow;
            outgoing[b] = outgoing[b] - flow;
            let resistor = match load {
                HarmonicElectrodeLoad::Resistive(r) => Some(&r.conductor),
                HarmonicElectrodeLoad::Reactive(ReactiveElectrodeLoad::SeriesRl {
                    conductor,
                    ..
                }) => conductor.as_ref(),
                HarmonicElectrodeLoad::Reactive(ReactiveElectrodeLoad::Capacitor { .. }) => None,
            };
            let average_heat_w = match resistor {
                Some(r) => {
                    0.5 * r
                        .joule_power_w(flow.abs())
                        .map_err(|_| invalid("harmonic conductor heat is not representable"))?
                }
                None => 0.0,
            };
            if !finite(flow)
                || !average_heat_w.is_finite()
                || (power(voltage, flow) - average_heat_w).abs() > controls.power_tolerance_w
            {
                return Err(invalid(
                    "harmonic load work does not match its source-owned heat",
                ));
            }
            heat += average_heat_w;
            load_responses.push(HarmonicElectrodeLoadResponse {
                load: load.clone(),
                voltage_v: voltage,
                current_a: flow,
                average_heat_w,
            });
        }
        let effective_conditions = source_conditions
            .iter()
            .zip(&outgoing)
            .map(|(&(v, z), &i)| (v, if v { z } else { z - charge(i, omega) }))
            .collect::<Vec<_>>();
        let mut parts = Vec::with_capacity(2);
        for (part, map) in boundary.iter().enumerate() {
            poll(cx)?;
            let scalar = |z: C64| if part == 0 { z.re } else { z.im };
            let full = full.iter().map(|&z| scalar(z)).collect::<Vec<_>>();
            let u = full[..nu]
                .chunks_exact(3)
                .map(|v| [v[0], v[1], v[2]])
                .collect();
            let voltage = full[nu..]
                .iter()
                .zip(&map.voltage_reference)
                .map(|(a, b)| a + b)
                .collect();
            let electrodes = electrodes
                .iter()
                .zip(&effective_conditions)
                .map(|(e, &(v, z))| PiezoelectricElectrode {
                    nodes: e.nodes,
                    condition: if v {
                        ElectrodeCondition::Voltage(scalar(z))
                    } else {
                        ElectrodeCondition::Charge(scalar(z))
                    },
                })
                .collect::<Vec<_>>();
            let f = effective_force
                .iter()
                .map(|f| f.map(scalar))
                .collect::<Vec<_>>();
            let q = reversible_charges
                .iter()
                .map(|&q| scalar(q))
                .collect::<Vec<_>>();
            parts.push(self.recover(
                &full,
                u,
                voltage,
                map,
                &electrodes,
                &f,
                &q,
                controls.field,
                iterations,
                relative_residual,
                cx,
            )?);
        }
        let (a, b) = (&parts[0], &parts[1]);
        let mut cells = Vec::with_capacity(self.tetrahedra.len());
        let mut material_heat = 0.0;
        let mut memory_energy = 0.0;
        for (id, (a, b)) in a.cells.iter().zip(&b.cells).enumerate() {
            poll(cx)?;
            let mut cell = HarmonicPiezoelectricCell {
                strain_mandel: combine(a.strain_mandel, b.strain_mandel),
                stress_mandel_pa: combine(a.response.stress_mandel_pa, b.response.stress_mandel_pa),
                electric_field_v_m: combine(
                    a.response.electric_field_v_m,
                    b.response.electric_field_v_m,
                ),
                displacement_c_m2: combine(
                    a.response.displacement_c_m2,
                    b.response.displacement_c_m2,
                ),
                average_mechanical_heat_w: 0.0,
                average_dielectric_heat_w: 0.0,
                average_memory_energy_j: 0.0,
            };
            if let Some(intrinsic) = &intrinsic {
                let response = &intrinsic.cells[id];
                let volume = self.dielectric.volumes_m3()[id];
                for axis in 0..6 {
                    cell.stress_mandel_pa[axis] =
                        cell.stress_mandel_pa[axis] + response.additional_stress_mandel_pa[axis];
                }
                for axis in 0..3 {
                    cell.displacement_c_m2[axis] =
                        cell.displacement_c_m2[axis] + response.polarization_c_m2[axis];
                }
                cell.average_mechanical_heat_w = volume * response.average_mechanical_heat_w_m3;
                cell.average_dielectric_heat_w = volume * response.average_dielectric_heat_w_m3;
                cell.average_memory_energy_j = volume
                    * (response.average_mechanical_memory_energy_j_m3
                        + response.average_dielectric_memory_energy_j_m3);
                material_heat += cell.average_mechanical_heat_w + cell.average_dielectric_heat_w;
                memory_energy += cell.average_memory_energy_j;
            }
            if cycle_norm(&cell.strain_mandel) > controls.field.maximum_strain_norm
                || cycle_norm(&cell.electric_field_v_m) > controls.field.maximum_electric_field_v_m
            {
                return Err(invalid(
                    "harmonic full-cycle field leaves the small-signal envelope",
                ));
            }
            cells.push(cell);
        }
        let displacement_m = a
            .displacement_m
            .iter()
            .zip(&b.displacement_m)
            .map(|(&a, &b)| combine(a, b))
            .collect::<Vec<_>>();
        let velocity_m_s = displacement_m
            .iter()
            .map(|u| u.map(|z| current(z, omega)))
            .collect::<Vec<_>>();
        let support_force_n = a
            .support_force_n
            .iter()
            .zip(&b.support_force_n)
            .map(|(&a, &b)| combine(a, b))
            .collect::<Vec<_>>();
        // Each real recovery checked one quadrature. Admission of a peak
        // phasor also needs its complex magnitude: independently allowing two
        // components up to the limit would weaken that limit by sqrt(2).
        let mut maximum_force_residual_n = 0.0_f64;
        let mut electric_residual = vec![C64::ZERO; boundary[0].rhs.len()];
        for node in 0..n {
            poll(cx)?;
            for axis in 0..3 {
                if boundary[0].index[3 * node + axis].is_some() {
                    maximum_force_residual_n =
                        maximum_force_residual_n.max(support_force_n[node][axis].abs());
                }
            }
            if let Some(index) = boundary[0].index[nu + node] {
                electric_residual[index] = electric_residual[index]
                    + C64::new(a.boundary_charge_c[node], b.boundary_charge_c[node]);
            }
        }
        for (electrode, &(voltage_controlled, prescribed_charge)) in
            electrodes.iter().zip(&effective_conditions)
        {
            poll(cx)?;
            if !voltage_controlled {
                let index = boundary[0].index[nu + electrode.nodes[0]]
                    .expect("admitted current electrode has one charge constraint");
                electric_residual[index] = electric_residual[index] - prescribed_charge;
            }
        }
        let mut maximum_charge_residual_c = 0.0_f64;
        for residual in electric_residual {
            poll(cx)?;
            if !finite(residual) {
                return Err(invalid("harmonic charge residual is not representable"));
            }
            maximum_charge_residual_c = maximum_charge_residual_c.max(residual.abs());
        }
        if maximum_force_residual_n > controls.field.force_tolerance_n {
            return Err(super::TetPiezoelectricError::Balance(
                "peak complex unconstrained force residual exceeds tolerance",
            ));
        }
        if maximum_charge_residual_c > controls.field.charge_tolerance_c {
            return Err(super::TetPiezoelectricError::Balance(
                "peak complex free-node or electrode-total charge residual exceeds tolerance",
            ));
        }
        let voltage_v = a
            .voltage_v
            .iter()
            .zip(&b.voltage_v)
            .map(|(&a, &b)| C64::new(a, b))
            .collect();
        let mut mechanical = 0.0;
        let mut electrical = 0.0;
        for node in 0..n {
            poll(cx)?;
            for axis in 0..3 {
                let reaction = if boundary[0].index[3 * node + axis].is_none() {
                    support_force_n[node][axis]
                } else {
                    C64::ZERO
                };
                mechanical += power(forces[node][axis] + reaction, velocity_m_s[node][axis]);
            }
            electrical += power(full[nu + node], current(charges[node], omega));
        }
        let mut electrode_responses = Vec::with_capacity(electrodes.len());
        for (id, e) in electrodes.iter().enumerate() {
            poll(cx)?;
            let charge_c = C64::new(a.electrodes[id].charge_c, b.electrodes[id].charge_c);
            let source_current_a = match e.condition {
                HarmonicElectrodeCondition::Current(i) => i,
                HarmonicElectrodeCondition::Voltage(_) => current(charge_c, omega) + outgoing[id],
            };
            electrical += power(full[nu + e.nodes[0]], source_current_a);
            electrode_responses.push(HarmonicElectrodeResponse {
                voltage_v: C64::new(a.electrodes[id].voltage_v, b.electrodes[id].voltage_v),
                charge_c,
                source_current_a,
            });
        }
        let storage = 0.5 * a.stored_energy_j + 0.5 * b.stored_energy_j + memory_energy;
        let residual = mechanical + electrical - heat - material_heat;
        if [
            kinetic,
            storage,
            mechanical,
            electrical,
            heat,
            material_heat,
            memory_energy,
            residual,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || velocity_m_s.iter().flatten().any(|&z| !finite(z))
        {
            return Err(invalid(
                "harmonic physical work/storage is not representable",
            ));
        }
        if residual.abs() > controls.power_tolerance_w {
            return Err(super::TetPiezoelectricError::Balance(
                "harmonic source power does not equal external load and intrinsic material heat",
            ));
        }
        poll(cx)?;
        Ok(HarmonicPiezoelectricResponse {
            displacement_m,
            velocity_m_s,
            voltage_v,
            support_force_n,
            cells,
            electrodes: electrode_responses,
            loads: load_responses,
            average_stored_energy_j: storage,
            average_kinetic_energy_j: kinetic,
            average_mechanical_power_w: mechanical,
            average_electrical_power_w: electrical,
            average_load_heat_w: heat,
            average_material_heat_w: material_heat,
            power_balance_residual_w: residual,
            maximum_force_residual_n,
            maximum_charge_residual_c,
            iterations,
            relative_residual,
        })
    }
}
