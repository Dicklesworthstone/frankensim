//! Retain several modes of each ORIGINAL regional material plate in one slit.
//! The existing plate/force-port reduction owns the eigenbasis, modal masses,
//! signed pressure areas, original force footprint and selected material loss.
use super::{AcousticRealizeError, ApertureMode, BladeState, ModalAperture,
    ModalApertureSpec, PairedApertureLimits, MAX_MODES, invalid};
use crate::thin_plate::{PlateChartRadiation, ResolvedPlateChart};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use crate::plate_viscosity::{PlateViscosity, PlateViscosityLimits};

/// Original plate reduction plus explicit spatial and work limits.
#[derive(Clone, Debug)]
pub struct ModalPlateOptions {
    /// Original eigenwindow, count, support, force footprint and loss selection.
    /// Selected material viscosity requires zero authored damping, as before.
    pub reduction: PlateChartRadiation,
    /// Unique physical boundary edges whose length-weighted mean opens the slit.
    pub slit_edges: Vec<[usize; 2]>,
    /// Maximum P1 displacement-trace slope of the SUMMED modal motion.
    /// This does not bound unexported DKT rotations or continuum derivatives.
    pub maximum_trace_slope: f64,
    /// Absolute departure of slit-node motion from its mean [m], per blade.
    /// This bounds the effective-uniform-slit approximation, not spatial contact.
    pub maximum_slit_deviation_m: f64,
    /// Maximum omega*dt for EVERY retained mode, in (0,1].
    pub maximum_angular_step: f64,
    /// Cold node ceiling, at most 512.
    pub max_nodes: usize,
    /// Cold triangle ceiling, at most 2048.
    pub max_triangles: usize,
    /// Original force-port projection work budget, including nodal shape ports.
    pub max_projection_terms: usize,
}

/// Source and physical displacement basis in original plate-node order.
#[derive(Clone, Debug)]
pub struct MaterialModalPlate {
    source: ResolvedPlateChart,
    options: ModalPlateOptions,
    first: usize,
    modes: Vec<ApertureMode>,
    shapes: Vec<Vec<f64>>,
    width_m: f64,
    slit_weights: Vec<f64>,
    viscosity: Option<PlateViscosity>,
    viscosity_matrix: Option<Vec<f64>>,
}
fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() { Err(AcousticRealizeError::Cancelled) } else { Ok(()) }
}
/// Admit the material adapter's declared angular cap before geometry work.
/// Actual modes also need the renderer's stricter phase-resolution admission.
pub fn admit_modal_angular_step(maximum_angular_step: f64, dt_s: f64)
    -> Result<(), AcousticRealizeError>
{
    if !maximum_angular_step.is_finite() || maximum_angular_step <= 0.0
        || maximum_angular_step > 1.0 || !dt_s.is_finite() || dt_s <= 0.0 {
        return Err(invalid("material modal angular cap must be in (0,1] radians with a positive finite timestep"));
    }
    Ok(())
}
impl MaterialModalPlate {
    /// Exact immutable regional material/geometry/loss source.
    pub const fn source(&self) -> &ResolvedPlateChart { &self.source }
    /// Original reduction, physical port and approximation declarations.
    pub const fn options(&self) -> &ModalPlateOptions { &self.options }
    /// Half-open coordinate range in the retained ModalAperture.
    pub fn mode_range(&self) -> std::ops::Range<usize> { self.first..self.first+self.modes.len() }
    /// One dimensionless displacement row per original node, in local mode order.
    /// These share the signed pressure/force basis; they are not separately solved.
    pub fn displacement_shapes(&self) -> &[Vec<f64>] { &self.shapes }
    /// Original length-integrated slit width [m]. It is not multiplied by modes.
    pub const fn width_m(&self) -> f64 { self.width_m }

    /// Original regional tensor selection; absent on the diagonal-loss image.
    pub fn viscosity(&self) -> Option<&PlateViscosity> { self.viscosity.as_ref() }
    /// Full local C [N s/m] in this plate's retained physical modal coordinates.
    pub fn viscosity_matrix_n_s_m(&self) -> Option<&[f64]> { self.viscosity_matrix.as_deref() }

    fn reduce(source: ResolvedPlateChart, options: ModalPlateOptions, first: usize,
        blade: usize, dt: f64, gate: &CancelGate,
        tensor: Option<(PlateViscosity, PlateViscosityLimits)>) -> Result<Self, AcousticRealizeError>
    {
        poll(gate)?;
        admit_modal_angular_step(options.maximum_angular_step, dt)?;
        let mesh = &source.chart().mesh;
        let n = mesh.nodes.len();
        if !(1..=512).contains(&options.max_nodes) || !(1..=2048).contains(&options.max_triangles)
            || n == 0 || n > options.max_nodes || mesh.tris.len() > options.max_triangles
            || options.reduction.n_modes == 0 || options.reduction.n_modes > MAX_MODES-first
            || options.slit_edges.is_empty() || options.slit_edges.len() > n
            || !options.maximum_trace_slope.is_finite() || options.maximum_trace_slope <= 0.0
            || !options.maximum_slit_deviation_m.is_finite() || options.maximum_slit_deviation_m < 0.0 {
            return Err(invalid("material modal plate requires bounded geometry, modes, clock and spatial limits"));
        }
        let terms = n.checked_mul(options.reduction.n_modes+1)
            .and_then(|v| mesh.tris.len().checked_mul(3).and_then(|t|v.checked_add(t)))
            .and_then(|v| v.checked_mul(n))
            .ok_or_else(||invalid("modal plate nodal projection work overflow"))?;
        if terms > options.max_projection_terms {
            return Err(invalid("modal plate nodal basis exceeds the original projection budget"));
        }
        let boundary: std::collections::BTreeSet<_> = mesh.boundary_edges().into_iter()
            .map(|(a,b)|(a.min(b),a.max(b))).collect();
        let mut seen = std::collections::BTreeSet::new();
        let mut weights = vec![0.0; n];
        let mut width = 0.0;
        for &[a,b] in &options.slit_edges {
            poll(gate)?;
            let edge = (a.min(b),a.max(b));
            if a >= n || b >= n || a == b || !boundary.contains(&edge) || !seen.insert(edge) {
                return Err(invalid("modal plate slit requires unique original boundary edges"));
            }
            let (a_pos,b_pos) = (mesh.nodes[a],mesh.nodes[b]);
            let length = (a_pos.0-b_pos.0).hypot(a_pos.1-b_pos.1);
            if !length.is_finite() || length <= 0.0 { return Err(invalid("modal plate slit length is invalid")); }
            width += length; weights[a] += 0.5*length; weights[b] += 0.5*length;
        }
        if !width.is_finite() || width <= 0.0 { return Err(invalid("modal slit width overflow")); }
        for w in &mut weights { *w /= width; }
        // Node unit-force columns reconstruct the ORIGINAL P1 displacement trace
        // from the very same eigensolve. Supported nodes remain exactly zero.
        let mut nodal_ports = vec![vec![0.0; n]; n];
        for (i,row) in nodal_ports.iter_mut().enumerate() { row[i] = 1.0; }
        poll(gate)?;
        let (bodies,ports,projected) = if let Some((law,limits)) = &tensor {
            let (image,ports) = source.certified_viscous_radiators_with_force_ports(
                &options.reduction, &nodal_ports, options.max_projection_terms, law, *limits, gate)?;
            (image.radiators().to_vec(), ports, Some(image.damping().matrix_per_s().to_vec()))
        } else {
            let (bodies,ports) = source.certified_radiators_with_force_ports(
                &options.reduction, &nodal_ports, options.max_projection_terms)?;
            (bodies,ports,None)
        };
        poll(gate)?;
        if bodies.len() != options.reduction.n_modes || ports.len() != n
            || ports.iter().any(|p|p.len() != bodies.len()) {
            return Err(invalid("requested material plate modes or physical shape ports were not resolved"));
        }
        let mut shapes = vec![vec![0.0; bodies.len()]; n];
        let mut modes = Vec::with_capacity(bodies.len());
        let mut coordinate_scale = Vec::with_capacity(bodies.len());
        for (i,body) in bodies.iter().enumerate() {
            poll(gate)?;
            let root_mass = body.mass_kg.sqrt();
            let scale = ports.iter().map(|p|(p[i]*root_mass).abs()).fold(0.0_f64,f64::max);
            if !scale.is_finite() || scale <= 0.0 || !body.omega.is_finite() || body.omega <= 0.0
                || body.omega*dt > options.maximum_angular_step {
                return Err(invalid("material mode has no finite displacement trace or exceeds its temporal resolution"));
            }
            for (row,port) in shapes.iter_mut().zip(&ports) { row[i] = port[i]*root_mass/scale; }
            let mass = body.mass_kg/scale/scale;
            coordinate_scale.push(root_mass/scale);
            let gap = weights.iter().zip(&shapes).map(|(w,row)|w*row[i]).sum();
            let mut force = [0.0;2]; force[blade] = body.drive_participation/scale;
            modes.push(ApertureMode { mass_kg:mass, stiffness_n_m:mass*body.omega*body.omega,
                damping_n_s_m:2.0*body.zeta*body.omega*mass,
                pressure_area_m2:body.area_m2/scale, gap_projection:gap, force_projection:force });
        }
        // Transform BOTH axes of C from the owner's mass-normalized coordinates.
        // The pressure/gap/force columns above came from that identical eigensolve.
        // Never divide by the mean slit displacement or discard signed cross terms.
        let viscosity_matrix = projected.map(|c| {
            let count = modes.len();
            let mut physical = vec![0.0; count*count];
            for i in 0..count { for j in i..count {
                let value = c[count*i+j]*coordinate_scale[i]*coordinate_scale[j];
                physical[count*i+j] = value; physical[count*j+i] = value;
            } }
            physical
        });
        Ok(Self { source, options, first, modes, shapes, width_m:width, slit_weights:weights,
            viscosity: tensor.map(|(law,_)|law), viscosity_matrix })
    }

    pub(super) fn validate_state(&self, all: &[BladeState]) -> Result<(), AcousticRealizeError> {
        let state = &all[self.mode_range()];
        let motion = |node:usize| self.shapes[node].iter().zip(state)
            .map(|(p,s)|p*s.displacement_m).sum::<f64>();
        let mean: f64 = self.modes.iter().zip(state).map(|(m,s)|m.gap_projection*s.displacement_m).sum();
        for (node,&w) in self.slit_weights.iter().enumerate() {
            if w > 0.0 && (motion(node)-mean).abs() > self.options.maximum_slit_deviation_m {
                return Err(invalid("summed material blade motion exceeds the effective-slit deviation limit"));
            }
        }
        let mesh = &self.source.chart().mesh;
        for tri in &mesh.tris {
            let [a,b,c] = tri.map(|i|mesh.nodes[i]);
            let [u,v,w] = tri.map(motion);
            let det = (b.0-a.0)*(c.1-a.1)-(b.1-a.1)*(c.0-a.0);
            let dx = ((v-u)*(c.1-a.1)-(w-u)*(b.1-a.1))/det;
            let dy = ((b.0-a.0)*(w-u)-(c.0-a.0)*(v-u))/det;
            if !dx.is_finite() || !dy.is_finite() || dx.hypot(dy) > self.options.maximum_trace_slope {
                return Err(invalid("summed material blade displacement exceeds its P1 trace-slope limit"));
            }
        }
        Ok(())
    }
}

impl ModalAperture {
    /// Derive multiple modes from one or two original material plates. Both
    /// plate frames point toward increasing total opening. Pressure, swept volume,
    /// mean gap, contact and each blade's original force footprint share one basis.
    /// No mode needs a nonzero mean gap or positive pressure projection by itself.
    ///
    /// Existing regional loss is retained by the shared diagonal weak-loss image;
    /// no hereditary or full cross-modal damping law is synthesized. There is one
    /// effective uniform slit/contact, NOT a distributed opposing-surface solve.
    /// P1 trace slope and slit deviation are checked on the SUMMED endpoint motion.
    /// Unexported DKT rotations and omitted higher modes are not certified here.
    /// Different source widths refuse; rest_opening_m is the TOTAL physical gap.
    /// Initial states address all requested modes, plate order then mode order.
    /// Cold cancellation brackets the existing non-interruptible modal reduction.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_plates(plates: Vec<(ResolvedPlateChart, ModalPlateOptions)>,
        rest_opening_m:f64, gas:GasState, impedance_pa_s_m3:f64, time_step_s:f64,
        max_steps:u64, limits:PairedApertureLimits, state:Vec<BladeState>,
        contact:Obstacle, gate:&CancelGate) -> Result<Self,AcousticRealizeError>
    {
        Self::from_material_images(plates.into_iter().map(|(s,o)|(s,o,None)).collect(),
            rest_opening_m,gas,impedance_pa_s_m3,time_step_s,max_steps,limits,state,contact,gate)
    }

    /// Select full source-derived anisotropic/regional viscosity for one or two
    /// material blades. Each tuple is (source, reduction, viscosity, loss limits).
    /// The existing plate owner projects C, force, pressure and every displacement
    /// trace on ONE basis per blade; the original coupled modal solver consumes
    /// the complete transformed C, not modal damping ratios or a waveform filter.
    ///
    /// Scalar material loss and authored damping must be absent. An explicit
    /// zero tensor can leave either blade lossless. Plates share fluid/contact
    /// coupling but do not acquire an invented solid-to-solid damping matrix.
    /// Source identity, tensor orientation and applicability stay attached.
    /// No viscosity is inferred from an elastic card or material name.
    ///
    /// All original geometry, clock, spatial, nonlinear and energy limits remain.
    /// This adds fixed-state Kelvin-Voigt loss, NOT Maxwell histories, distributed
    /// closure, curved-shell physics or a measured timbre claim.
    #[allow(clippy::too_many_arguments)]
    pub fn from_viscous_material_plates(
        plates: Vec<(ResolvedPlateChart, ModalPlateOptions, PlateViscosity, PlateViscosityLimits)>,
        rest_opening_m:f64, gas:GasState, impedance_pa_s_m3:f64, time_step_s:f64,
        max_steps:u64, limits:PairedApertureLimits, state:Vec<BladeState>,
        contact:Obstacle, gate:&CancelGate) -> Result<Self,AcousticRealizeError>
    {
        Self::from_material_images(plates.into_iter().map(|(s,o,c,l)|(s,o,Some((c,l)))).collect(),
            rest_opening_m,gas,impedance_pa_s_m3,time_step_s,max_steps,limits,state,contact,gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn from_material_images(
        plates: Vec<(ResolvedPlateChart, ModalPlateOptions, Option<(PlateViscosity, PlateViscosityLimits)>)>,
        rest_opening_m:f64, gas:GasState, impedance_pa_s_m3:f64, time_step_s:f64,
        max_steps:u64, limits:PairedApertureLimits, state:Vec<BladeState>,
        contact:Obstacle, gate:&CancelGate) -> Result<Self,AcousticRealizeError>
    {
        poll(gate)?;
        gas.validate_coherence().map_err(super::nonlinear)?;
        if !(1..=2).contains(&plates.len()) { return Err(invalid("modal valve needs one or two original plates")); }
        let count = plates.iter().try_fold(0usize,|n,(_,o,_)|n.checked_add(o.reduction.n_modes))
            .ok_or_else(||invalid("material modal count overflow"))?;
        if count == 0 || count > MAX_MODES || state.len() != count {
            return Err(invalid("material modal state must match every requested mode within the shared limit"));
        }
        let mut retained = Vec::with_capacity(plates.len());
        let mut modes = Vec::with_capacity(count);
        for (blade,(source,options,tensor)) in plates.into_iter().enumerate() {
            let plate = MaterialModalPlate::reduce(source,options,modes.len(),blade,time_step_s,gate,tensor)?;
            modes.extend_from_slice(&plate.modes); retained.push(plate);
        }
        let width = retained[0].width_m;
        if retained.iter().any(|p|p.width_m.to_bits() != width.to_bits()) {
            return Err(invalid("opposing modal plates require the same physical slit width"));
        }
        let mut valve = Self::new(ModalApertureSpec {rest_opening_m,width_m:width,modes,
            density_kg_m3:gas.density,impedance_pa_s_m3,time_step_s,max_steps,limits},state,contact)?;
        if retained.iter().any(|p|p.viscosity_matrix.is_some()) {
            let mut matrix = vec![0.0; count*count];
            for plate in &retained {
                let c = plate.viscosity_matrix.as_ref().ok_or_else(||
                    invalid("full modal loss requires an explicit tensor selection on every blade"))?;
                let n = plate.modes.len();
                for i in 0..n { for j in 0..n {
                    matrix[count*(plate.first+i)+plate.first+j] = c[n*i+j];
                } }
            }
            valve = valve.with_viscosity(matrix,gate)?;
        }
        valve.material_plates = retained; valve.material_gas = Some(gas);
        valve.validate_state(&valve.state)?;
        poll(gate)?;
        Ok(valve)
    }
    /// Original material plates, shapes and projection ranges, not copied state.
    pub fn material_plates(&self) -> &[MaterialModalPlate] { &self.material_plates }
    /// Exact material gas on the plate-driven path; absent on the raw modal path.
    pub const fn gas(&self) -> Option<&GasState> { self.material_gas.as_ref() }
    /// Actual nodal transverse displacement and velocity [m,m/s] on one source.
    /// This is a read of the retained signed basis, never another vibration solve.
    pub fn nodal_motion(&self,plate:usize,node:usize) -> Result<BladeState,AcousticRealizeError> {
        let plate = self.material_plates.get(plate).ok_or_else(||invalid("no material plate at this address"))?;
        let shape = plate.shapes.get(node).ok_or_else(||invalid("no material node at this address"))?;
        let state = &self.state[plate.mode_range()];
        let result = BladeState {
            displacement_m:shape.iter().zip(state).map(|(p,s)|p*s.displacement_m).sum(),
            velocity_m_s:shape.iter().zip(state).map(|(p,s)|p*s.velocity_m_s).sum(),
        };
        if !result.displacement_m.is_finite() || !result.velocity_m_s.is_finite() {
            return Err(invalid("material node motion overflow"));
        }
        Ok(result)
    }
}
