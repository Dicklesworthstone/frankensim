//! Tensioned filaments colliding with a sampled moving surface.
//!
//! Geometry -> the existing prestressed-beam frequencies and mass-normalized
//! sine shapes -> an ImpactBody and an fs-dcontact obstacle. No new time
//! integrator, synthesized buzz, random excitation, or sound-pressure gain.
//! Contact is integrated along the span: quadrature weights are METRES, not
//! 1/N. Stiffness is per unit span and must not change on contact refinement.
//!
//! A coiled strand may be homogenized with explicitly supplied effective
//! tension/bending data. Helix geometry below identifies only its mass, NOT
//! those constitutive properties. Small transverse motion, fixed endpoints,
//! one polarization; no coil torsion, end-plate motion or inter-wire friction.

use super::{BodyPotential, ImpactBody, ImpactError};
use crate::modal_acoustic_time::{ModalAcousticState, MAX_TIME_DOMAIN_ACOUSTIC_MODES};
use fs_dcontact::{Obstacle, string_collocation};
use fs_plate::{ModePair, shell::head::TensionedDisk};
use std::ops::Range;

fn invalid(what: &'static str) -> ImpactError { ImpactError::Invalid(what) }

/// Geometry of a uniform helical wire's centreline and circular metal section.
#[derive(Clone, Copy, Debug)]
pub struct HelicalWire {
    /// Radius of the metal wire [m], not the outside radius of the coil.
    pub wire_radius_m: f64,
    /// Helix centreline radius [m]; zero is a straight wire.
    pub coil_radius_m: f64,
    /// Axial advance per complete turn [m].
    pub pitch_m: f64,
    /// Density of the actual wire material [kg/m^3].
    pub density_kg_m3: f64,
}
impl HelicalWire {
    /// Material mass per AXIAL metre, including the helix's extra arc length.
    /// No tension, flexural stiffness or dissipative coefficient is inferred.
    pub fn linear_density_kg_m(self) -> Result<f64, ImpactError> {
        if [self.wire_radius_m, self.pitch_m, self.density_kg_m3].iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || !self.coil_radius_m.is_finite() || self.coil_radius_m < 0.0
        { return Err(invalid("wire section/density/pitch must be positive finite; coil radius nonnegative")); }
        let arc_ratio = (std::f64::consts::TAU * self.coil_radius_m / self.pitch_m).hypot(1.0);
        let density = self.density_kg_m3 * std::f64::consts::PI
            * self.wire_radius_m * self.wire_radius_m * arc_ratio;
        if !density.is_finite() || density <= 0.0 { return Err(invalid("helical line density is unrepresentable")); }
        Ok(density)
    }
}

/// Uniform fixed-end transverse filament in the receiver's reference XY plane.
/// The endpoint distance is its speaking length; both coordinates move in the
/// same declared transverse direction when coupled by `contact`.
#[derive(Clone, Debug)]
pub struct WireSpan {
    /// Physical endpoints [m]. End-plate compliance must be modeled separately.
    pub endpoints_m: [[f64; 2]; 2],
    /// Mass per axial length [kg/m]. May come from a measured coil or HelicalWire.
    pub linear_density_kg_m: f64,
    /// Installed effective axial tension [N]; not a pitch parameter.
    pub tension_n: f64,
    /// Effective transverse flexural stiffness [N m^2]. Zero is an ideal string.
    pub bending_stiffness_n_m2: f64,
    /// Explicit viscous modal coefficients [1/s]; one per retained sine mode.
    pub damping_per_s: Vec<f64>,
}
impl WireSpan {
    /// Length in the declared reference plane [m].
    #[must_use]
    pub fn length_m(&self) -> f64 {
        (self.endpoints_m[1][0]-self.endpoints_m[0][0])
            .hypot(self.endpoints_m[1][1]-self.endpoints_m[0][1])
    }
    fn validate(&self) -> Result<(), ImpactError> {
        if self.endpoints_m.iter().flatten().any(|v| !v.is_finite())
            || [self.length_m(), self.linear_density_kg_m, self.tension_n]
                .iter().any(|v| !v.is_finite() || *v <= 0.0)
            || !self.bending_stiffness_n_m2.is_finite() || self.bending_stiffness_n_m2 < 0.0
            || self.damping_per_s.is_empty() || self.damping_per_s.len() > MAX_TIME_DOMAIN_ACOUSTIC_MODES
            || self.damping_per_s.iter().any(|d| !d.is_finite() || *d < 0.0)
        { return Err(invalid("wire needs finite positive length/mass/tension and bounded explicit modal damping")); }
        Ok(())
    }
    /// Build existing linear mechanics with caller-supplied initial motion.
    /// The runtime owner additionally checks its clock and state/energy limits.
    pub fn body(&self, initial: Vec<ModalAcousticState>) -> Result<ImpactBody, ImpactError> {
        self.validate()?;
        if initial.len() != self.damping_per_s.len() || initial.iter().any(|s|
            !s.displacement_m_sqrt_kg.is_finite() || !s.velocity_m_sqrt_kg_per_s.is_finite())
        { return Err(invalid("wire initial motion must match every retained mass-normalized mode")); }
        let omega: Vec<_> = (1..=initial.len()).map(|k| fs_nlmodal::prestressed_beam_omega(
            self.length_m(), self.tension_n, self.linear_density_kg_m,
            self.bending_stiffness_n_m2, k,
        )).collect();
        if omega.iter().any(|w| !w.is_finite() || *w <= 0.0) {
            return Err(invalid("wire derived frequencies are unrepresentable"));
        }
        Ok(ImpactBody { potential: BodyPotential::Linear(omega), initial,
            damping_per_s: self.damping_per_s.clone() })
    }
    /// Fixed-order physical XY positions at the contact's declared stations.
    pub fn positions(&self, line: &LineContact) -> Result<Vec<[f64; 2]>, ImpactError> {
        self.validate()?;
        line.validate(self.length_m())?;
        Ok(line.stations_m.iter().map(|s| {
            let t = s/self.length_m();
            std::array::from_fn(|j| (1.0-t)*self.endpoints_m[0][j] + t*self.endpoints_m[1][j])
        }).collect())
    }
    /// Assemble closure = receiver displacement - wire displacement - gap.
    /// The SAME signed rows spread equal/opposite contact work to both bodies.
    /// `receiver_shapes` is point-major and may come from any geometric solver.
    /// Neither the initial state nor an already admitted model is mutated.
    pub fn contact(&self, line: &LineContact, receiver_shapes: &[Vec<f64>],
        receiver_modes: Range<usize>, wire_modes: Range<usize>, total_modes: usize,
    ) -> Result<Obstacle, ImpactError> {
        self.validate()?;
        line.validate(self.length_m())?;
        if total_modes == 0 || total_modes > MAX_TIME_DOMAIN_ACOUSTIC_MODES
            || receiver_modes.start >= receiver_modes.end || receiver_modes.end > total_modes
            || wire_modes.start >= wire_modes.end || wire_modes.end > total_modes
            || wire_modes.len() != self.damping_per_s.len()
            || (receiver_modes.start < wire_modes.end && wire_modes.start < receiver_modes.end)
            || receiver_shapes.len() != line.stations_m.len()
            || receiver_shapes.iter().any(|r| r.len() != receiver_modes.len() || r.iter().any(|v| !v.is_finite()))
        { return Err(invalid("wire contact requires disjoint complete mode ranges and finite point-major receiver shapes")); }
        let wire_shapes = string_collocation(self.length_m(), self.linear_density_kg_m,
            &line.stations_m, wire_modes.len()).map_err(|e| ImpactError::Owner(e.to_string()))?;
        let mut columns = vec![0.0; total_modes*line.stations_m.len()];
        for (i, row) in columns.chunks_exact_mut(total_modes).enumerate() {
            row[receiver_modes.clone()].copy_from_slice(&receiver_shapes[i]);
            for (j, slot) in row[wire_modes.clone()].iter_mut().enumerate() {
                *slot = -wire_shapes[i*wire_modes.len()+j];
            }
        }
        Obstacle::new(columns, line.stations_m.len(), total_modes, line.gaps_m.clone(),
            line.measures_m.clone(), line.stiffness_per_length, line.alpha, line.provenance.clone())
            .and_then(|o| o.with_internal_loss(line.internal_loss_s_m))
            .map_err(|e| ImpactError::Owner(e.to_string()))
    }
}

/// A line-distributed constitutive law with explicit spatial quadrature.
#[derive(Clone, Debug)]
pub struct LineContact {
    /// Increasing positions along the speaking length, strictly inside (0,L) [m].
    pub stations_m: Vec<f64>,
    /// Positive quadrature lengths [m]. They are not normalized to one.
    pub measures_m: Vec<f64>,
    /// Rest clearance [m] at each station; negative means declared interference.
    pub gaps_m: Vec<f64>,
    /// K in force-per-length = K [penetration]^alpha [N/m^(alpha+1)].
    pub stiffness_per_length: f64,
    /// Power-law exponent, at least one.
    pub alpha: f64,
    /// Existing nonadhesive Hunt-Crossley coefficient [s/m].
    pub internal_loss_s_m: f64,
    /// Physical/source or explicit estimate identity; never an authority claim.
    pub provenance: String,
}
impl LineContact {
    /// Uniform midpoint quadrature over the complete span, with supplied gap.
    /// This does not assert adequate contact resolution; refine it independently
    /// of the retained wire and surface modal bases.
    pub fn uniform(length_m: f64, cells: usize, gap_m: f64,
        stiffness_per_length: f64, alpha: f64, internal_loss_s_m: f64, provenance: String,
    ) -> Result<Self, ImpactError> {
        if !length_m.is_finite() || length_m <= 0.0 || !(1..=4096).contains(&cells) {
            return Err(invalid("line quadrature requires positive finite length and 1..=4096 cells"));
        }
        let dx = length_m/cells as f64;
        let line = Self { stations_m: (0..cells).map(|i| (i as f64+0.5)*dx).collect(),
            measures_m: vec![dx; cells], gaps_m: vec![gap_m; cells],
            stiffness_per_length, alpha, internal_loss_s_m, provenance };
        line.validate(length_m)?;
        Ok(line)
    }
    fn validate(&self, length: f64) -> Result<(), ImpactError> {
        let n = self.stations_m.len();
        if n == 0 || n > 4096 || self.measures_m.len() != n || self.gaps_m.len() != n
            || self.stations_m.iter().enumerate().any(|(i,s)| !s.is_finite() || *s <= 0.0
                || *s >= length || (i > 0 && *s <= self.stations_m[i-1]))
            || self.measures_m.iter().any(|w| !w.is_finite() || *w <= 0.0)
            || self.gaps_m.iter().any(|g| !g.is_finite())
            || !self.stiffness_per_length.is_finite() || self.stiffness_per_length < 0.0
            || !self.alpha.is_finite() || self.alpha < 1.0
            || !self.internal_loss_s_m.is_finite() || self.internal_loss_s_m < 0.0
            || self.provenance.trim().is_empty()
        { return Err(invalid("line contact requires finite ordered stations, positive lengths, gaps and an explicit law")); }
        Ok(())
    }
}

/// Sample the existing film's P1 transverse field at ACTUAL XY coordinates.
/// No nearest-node snapping: barycentric interpolation and its force transpose
/// use exactly the same point. Shared-edge ties choose the first mesh facet.
/// This is the existing linear transverse interpolation, not full DKT interior
/// recovery or a geometric/eigenbasis adequacy certificate.
pub fn film_shapes(film: &TensionedDisk, modes: &[ModePair], points: &[[f64; 2]])
    -> Result<Vec<Vec<f64>>, ImpactError>
{
    if modes.is_empty() || modes.len() > MAX_TIME_DOMAIN_ACOUSTIC_MODES
        || points.is_empty() || points.len() > 4096
        || modes.iter().any(|m| m.phi.len() != film.model.free || m.phi.iter().any(|v| !v.is_finite()))
        || points.iter().flatten().any(|p| !p.is_finite())
    { return Err(invalid("film contact samples require finite points and modes in the original pencil")); }
    let mut result = Vec::with_capacity(points.len());
    for &[x,y] in points {
        let mut found = None;
        for tri in &film.mesh.tris {
            let (ax,ay) = film.mesh.nodes[tri[0]];
            let (bx,by) = film.mesh.nodes[tri[1]];
            let (cx,cy) = film.mesh.nodes[tri[2]];
            let determinant = (bx-ax)*(cy-ay)-(by-ay)*(cx-ax);
            if !determinant.is_finite() || determinant <= 0.0 {
                return Err(invalid("film contact requires finite positively oriented facets"));
            }
            let b = ((x-ax)*(cy-ay)-(y-ay)*(cx-ax))/determinant;
            let c = ((bx-ax)*(y-ay)-(by-ay)*(x-ax))/determinant;
            let a = 1.0-b-c;
            // Only roundoff at shared edges is tolerated; there is no spatial
            // snapping radius and no extrapolation beyond the mesh polygon.
            if [a,b,c].iter().all(|v| v.is_finite() && *v >= -32.0*f64::EPSILON) {
                let mut bary = [a.max(0.0), b.max(0.0), c.max(0.0)];
                let sum = bary.iter().sum::<f64>();
                for v in &mut bary { *v /= sum; }
                found = Some((tri, bary));
                break;
            }
        }
        let (tri,bary) = found.ok_or_else(|| invalid("wire contact station lies outside the actual film mesh"))?;
        let row: Vec<_> = modes.iter().map(|mode| tri.iter().zip(bary).map(|(&node,b)|
            b*film.model.dof_map[3*node].map_or(0.0, |i| mode.phi[i])).sum::<f64>()).collect();
        if row.iter().any(|v| !v.is_finite()) { return Err(invalid("film contact interpolation overflow")); }
        result.push(row);
    }
    Ok(result)
}
