//! Reflected shock loading on stationary, explicitly matched structural surfaces.
//!
//! The original gas owner supplies its accepted wall impulse. The Euler path
//! retains pressure; the molecular paths include viscous normal stress and
//! resolved no-slip wall shear from both molecular half steps. Neither samples
//! endpoint pressure or replaces
//! the shock with a linear acoustic wave. The original
//! shell pressure port owns the spatial P1 load and its conjugate volume row.
//! Cold preparation compresses those rows; interval projection visits only
//! retained nodes. No gas or structural state is advanced by this adapter.
//!
//! Forces are the held interval means. Their impulse and spatial virtual work
//! are preserved; unresolved time variation still requires a timestep ladder.
//! A stationary wall does no gas work. Driving a deformable body with these
//! loads is one-way prescribed loading, not a reciprocal moving-wall FSI solve.

use fs_gas::dynamics::{Boundaries, Boundary, Euler1d, StepReport};
use fs_math::dd::Dd;
use fs_plate::shell::{
    ShellMesh,
    pressure::{OrientedShellFacet, ShellFace, ShellPressurePort},
    profile::ProfileBudget,
};
use std::collections::BTreeSet;

#[cfg(feature = "experimental-compressible")]
mod cartesian;

/// Input, geometry and cancellation refusals publish no partial load.
#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    /// Inconsistent geometry, boundary, clock or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The original shell pressure owner refused the surface.
    Surface(fs_plate::PlateError),
    /// Cancellation before returning a complete result.
    Cancelled,
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Euler wall/structural loading: {self:?}")
    }
}
impl std::error::Error for Error {}
impl From<fs_plate::PlateError> for Error {
    fn from(error: fs_plate::PlateError) -> Self {
        if error == fs_plate::PlateError::Cancelled {
            Self::Cancelled
        } else {
            Self::Surface(error)
        }
    }
}

/// End of a one-dimensional branch, in increasing axial-coordinate order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum End {
    /// Smallest axial coordinate.
    Left,
    /// Largest axial coordinate.
    Right,
}

/// Actual accepted pressure-time integral on one rigid fluid boundary.
/// Fields are private so a transported momentum flux cannot masquerade as
/// reflected-wall pressure. Exterior pressure is absolute, held over the SAME
/// accepted interval; it is never subtracted from the gas owner's flux.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallPressureImpulse {
    inner: WallImpulse,
}

/// Signed normal traction from a complete accepted molecular/Euler interval.
///
/// The gas contribution is the integral of `p - tau_nn`, not absolute pressure.
/// It may be tensile. Tangential slip-wall traction, wall work and heat must
/// remain zero; only the original combined molecular owner can construct it.
#[cfg(feature = "experimental-compressible")]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallTractionImpulse {
    inner: WallImpulse,
}

/// Complete vector traction from one accepted stationary molecular wall face.
///
/// The gas impulse includes pressure and all resolved viscous momentum transfer,
/// including shear. Only the original complete molecular flow owner can create
/// this account. Open, moving and nonadiabatic walls require additional mass or
/// energy coupling and are not admitted. The structural frame is supplied by an
/// explicit orthonormal embedding; the gas face is never spatially averaged.
#[cfg(feature = "experimental-compressible")]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallVectorTractionImpulse {
    normal: WallImpulse,
    tangential_pa_s: [f64; 3],
}

/// Common stationary normal-load account. Only the pressure wrapper imposes
/// positive gas pressure; molecular normal stress is signed without clipping.
#[derive(Clone, Copy, Debug, PartialEq)]
struct WallImpulse {
    duration_s: f64,
    end_time_s: f64,
    area_m2: f64,
    normal: [f64; 3],
    gas_pa_s: f64,
    exterior_pa_s: f64,
}
impl WallPressureImpulse {
    /// Consume the latest accepted frozen-gas report at a reflecting end.
    /// `axis` is the unit vector of increasing branch coordinate in the same
    /// global frame as the structural surface. Both SSP RK stages contribute.
    /// The caller supplies its latest report; current extensive state is
    /// checked, but these public reports are not authenticated history handles.
    pub fn from_euler1d(
        flow: &Euler1d,
        report: &StepReport,
        end: End,
        axis: [f64; 3],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let boundary = match (flow.boundaries(), end) {
            (Boundaries::Ends { left, .. }, End::Left) => left,
            (Boundaries::Ends { right, .. }, End::Right) => right,
            (Boundaries::Periodic, _) => {
                return Err(Error::Invalid("periodic faces have no wall load"));
            }
        };
        if boundary != Boundary::Reflecting || report.after != flow.totals() {
            return Err(Error::Invalid(
                "latest accepted reflecting-wall report required",
            ));
        }
        let (change, face, sign) = match end {
            End::Left => (report.left_boundary_change, 0, -1.0),
            End::Right => (report.right_boundary_change, flow.mesh().cell_count(), 1.0),
        };
        if change.mass != 0.0 || change.energy != 0.0 {
            return Err(Error::Invalid(
                "stationary adiabatic wall must transfer only momentum",
            ));
        }
        Self::new(
            report.elapsed.value(),
            flow.time().value(),
            flow.mesh().face_areas()[face],
            axis.map(|value| sign * value),
            -sign * change.momentum,
            exterior_pressure_pa,
        )
    }

    /// Same accepted wall account for every species in the original mixture.
    /// No constant-gamma substitution or separate species/energy transport is
    /// introduced. The composition-aware wall solver already owns reflection.
    pub fn from_mixture_euler1d(
        flow: &fs_gas::multicomponent::MixtureEuler1d,
        report: &fs_gas::multicomponent::StepReport,
        end: End,
        axis: [f64; 3],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let boundary = match (flow.boundaries(), end) {
            (Boundaries::Ends { left, .. }, End::Left) => left,
            (Boundaries::Ends { right, .. }, End::Right) => right,
            (Boundaries::Periodic, _) => {
                return Err(Error::Invalid("periodic faces have no wall load"));
            }
        };
        if boundary != Boundary::Reflecting || report.after != flow.totals() {
            return Err(Error::Invalid(
                "latest accepted reflecting-wall mixture report required",
            ));
        }
        let (change, face, sign) = match end {
            End::Left => (report.left_boundary_change, 0, -1.0),
            End::Right => (report.right_boundary_change, flow.mesh().cell_count(), 1.0),
        };
        if change.species_masses.iter().any(|mass| *mass != 0.0) || change.energy != 0.0 {
            return Err(Error::Invalid(
                "stationary adiabatic wall must transfer only momentum",
            ));
        }
        Self::new(
            report.elapsed.value(),
            flow.time().value(),
            flow.mesh().face_areas()[face],
            axis.map(|value| sign * value),
            -sign * change.momentum,
            exterior_pressure_pa,
        )
    }

    fn new(
        duration_s: f64,
        end_time_s: f64,
        area_m2: f64,
        normal: [f64; 3],
        gas_impulse_n_s: f64,
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        if !positive(gas_impulse_n_s) {
            return Err(Error::Invalid(
                "reflected gas pressure impulse must be finite and positive",
            ));
        }
        Ok(Self {
            inner: WallImpulse::new(
                duration_s,
                end_time_s,
                area_m2,
                normal,
                gas_impulse_n_s,
                exterior_pressure_pa,
            )?,
        })
    }

    /// Actual accepted duration [s].
    pub const fn duration_s(&self) -> f64 {
        self.inner.duration_s
    }
    /// Accepted gas clock at the end of the interval [s].
    pub const fn end_time_s(&self) -> f64 {
        self.inner.end_time_s
    }
    /// Physical gas face area [m²].
    pub const fn area_m2(&self) -> f64 {
        self.inner.area_m2
    }
    /// Positive unit normal from the gas into the support.
    pub const fn outward_normal(&self) -> [f64; 3] {
        self.inner.normal
    }
    /// Original wall pressure integrated over both accepted RK stages [Pa s].
    pub const fn gas_pressure_impulse_pa_s(&self) -> f64 {
        self.inner.gas_pa_s
    }
    /// Opposing held exterior absolute pressure times accepted duration [Pa s].
    pub const fn exterior_pressure_impulse_pa_s(&self) -> f64 {
        self.inner.exterior_pa_s
    }
    /// Signed gas-minus-exterior pressure-time integral [Pa s].
    pub fn net_pressure_impulse_pa_s(&self) -> f64 {
        self.inner.net_pa_s()
    }
    /// Actual interval-mean reflected gas pressure [Pa].
    pub fn mean_gas_pressure_pa(&self) -> f64 {
        self.inner.gas_pa_s / self.inner.duration_s
    }
}

#[cfg(feature = "experimental-compressible")]
impl WallTractionImpulse {
    fn new(
        duration_s: f64,
        end_time_s: f64,
        area_m2: f64,
        normal: [f64; 3],
        gas_impulse_n_s: f64,
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        Ok(Self {
            inner: WallImpulse::new(
                duration_s,
                end_time_s,
                area_m2,
                normal,
                gas_impulse_n_s,
                exterior_pressure_pa,
            )?,
        })
    }
    /// Complete accepted molecular/Euler duration [s].
    pub const fn duration_s(&self) -> f64 {
        self.inner.duration_s
    }
    /// Complete accepted source clock at the end of the interval [s].
    pub const fn end_time_s(&self) -> f64 {
        self.inner.end_time_s
    }
    /// Original physical gas face area [m²].
    pub const fn area_m2(&self) -> f64 {
        self.inner.area_m2
    }
    /// Unit normal from the gas into the support, even for tensile traction.
    pub const fn outward_normal(&self) -> [f64; 3] {
        self.inner.normal
    }
    /// Signed integral of gas `p - tau_nn` over the full split interval [Pa s].
    pub const fn gas_traction_impulse_pa_s(&self) -> f64 {
        self.inner.gas_pa_s
    }
    /// Held opposing exterior absolute pressure times accepted duration [Pa s].
    pub const fn exterior_pressure_impulse_pa_s(&self) -> f64 {
        self.inner.exterior_pa_s
    }
    /// Signed gas normal traction minus exterior pressure, integrated [Pa s].
    pub fn net_traction_impulse_pa_s(&self) -> f64 {
        self.inner.net_pa_s()
    }
    /// Interval-mean gas normal traction [Pa]; this is not absolute pressure.
    pub fn mean_gas_traction_pa(&self) -> f64 {
        self.inner.gas_pa_s / self.inner.duration_s
    }
}

#[cfg(feature = "experimental-compressible")]
impl WallVectorTractionImpulse {
    fn new(
        duration_s: f64,
        end_time_s: f64,
        area_m2: f64,
        normal: [f64; 3],
        normal_impulse_n_s: f64,
        tangential_impulse_n_s: [f64; 3],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let normal = WallImpulse::new(
            duration_s,
            end_time_s,
            area_m2,
            normal,
            normal_impulse_n_s,
            exterior_pressure_pa,
        )?;
        if tangential_impulse_n_s.iter().any(|x| !x.is_finite()) {
            return Err(Error::Invalid(
                "finite original tangential wall impulse required",
            ));
        }
        let tangential_pa_s = tangential_impulse_n_s.map(|x| x / area_m2);
        let result = Self {
            normal,
            tangential_pa_s,
        };
        if (0..3).any(|axis| {
            !tangential_pa_s[axis].is_finite()
                || (tangential_impulse_n_s[axis] != 0.0 && tangential_pa_s[axis] == 0.0)
        }) || result
            .gas_traction_impulse_pa_s()
            .iter()
            .any(|x| !x.is_finite())
            || result
                .net_traction_impulse_pa_s()
                .iter()
                .any(|x| !x.is_finite())
        {
            return Err(Error::Invalid(
                "wall vector-traction integral is unrepresentable",
            ));
        }
        Ok(result)
    }

    /// Complete accepted molecular/Euler duration [s].
    pub const fn duration_s(&self) -> f64 {
        self.normal.duration_s
    }
    /// Complete accepted source clock at the end of the interval [s].
    pub const fn end_time_s(&self) -> f64 {
        self.normal.end_time_s
    }
    /// Original physical gas face area [m²].
    pub const fn area_m2(&self) -> f64 {
        self.normal.area_m2
    }
    /// Unit normal from the gas into the support in the structural frame.
    pub const fn outward_normal(&self) -> [f64; 3] {
        self.normal.normal
    }
    /// Complete gas traction-time vector in the structural frame [Pa s].
    pub fn gas_traction_impulse_pa_s(&self) -> [f64; 3] {
        core::array::from_fn(|axis| {
            self.normal.normal[axis] * self.normal.gas_pa_s + self.tangential_pa_s[axis]
        })
    }
    /// Opposing exterior absolute pressure integrated over this interval [Pa s].
    pub const fn exterior_pressure_impulse_pa_s(&self) -> f64 {
        self.normal.exterior_pa_s
    }
    /// Gas vector traction minus the exterior normal pressure, integrated [Pa s].
    pub fn net_traction_impulse_pa_s(&self) -> [f64; 3] {
        // Subtract opposing normal pressures before combining them with
        // shear. A large common ambient pressure must not erase a small drag.
        let normal = self.normal.net_pa_s();
        core::array::from_fn(|axis| self.normal.normal[axis] * normal + self.tangential_pa_s[axis])
    }
    /// Interval-mean gas traction vector [Pa], including resolved wall shear.
    pub fn mean_gas_traction_pa(&self) -> [f64; 3] {
        self.gas_traction_impulse_pa_s()
            .map(|x| x / self.normal.duration_s)
    }
}

impl WallImpulse {
    fn new(
        duration_s: f64,
        end_time_s: f64,
        area_m2: f64,
        normal: [f64; 3],
        gas_impulse_n_s: f64,
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        if !positive(duration_s)
            || !end_time_s.is_finite()
            || end_time_s < duration_s
            || !positive(area_m2)
            || !gas_impulse_n_s.is_finite()
            || !unit(normal)
            || !exterior_pressure_pa.is_finite()
            || exterior_pressure_pa < 0.0
        {
            return Err(Error::Invalid(
                "finite wall interval, area, normal load and unit frame required",
            ));
        }
        let gas_pa_s = gas_impulse_n_s / area_m2;
        let exterior_pa_s = exterior_pressure_pa * duration_s;
        if !gas_pa_s.is_finite()
            || (gas_impulse_n_s != 0.0 && gas_pa_s == 0.0)
            || !exterior_pa_s.is_finite()
            || (exterior_pressure_pa > 0.0 && exterior_pa_s == 0.0)
            || !(gas_pa_s - exterior_pa_s).is_finite()
        {
            return Err(Error::Invalid(
                "wall normal-traction integral is unrepresentable",
            ));
        }
        Ok(Self {
            duration_s,
            end_time_s,
            area_m2,
            normal,
            gas_pa_s,
            exterior_pa_s,
        })
    }

    fn net_pa_s(self) -> f64 {
        self.gas_pa_s - self.exterior_pa_s
    }
}

/// Cold preparation limits include the original global indexing and all dense
/// pressure rows visited before their sparse compression. No mesh is truncated.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceBudget {
    /// Original shell geometry limits, passed unchanged to its pressure owner.
    pub geometry: ProfileBudget,
    /// Maximum separately loaded fluid faces.
    pub maximum_patches: usize,
    /// Maximum dense-row and selected-geometry scalar visits during preparation.
    pub maximum_preparation_terms: usize,
}

#[derive(Clone, Debug)]
struct PatchWeight {
    node: usize,
    normal_area_m2: [f64; 3],
    #[cfg(feature = "experimental-compressible")]
    area_m2: f64,
}

#[derive(Clone, Debug)]
struct Patch {
    area_m2: f64,
    normal: [f64; 3],
    row: Vec<PatchWeight>,
}

/// One planar oriented structural patch per fluid face, preserving original
/// node indices. Selected triangles may not be repeated across fluid faces.
#[derive(Clone, Debug)]
pub struct ShockLoadSurface {
    positions: Vec<[f64; 3]>,
    patches: Vec<Patch>,
}
impl ShockLoadSurface {
    /// Prepare the existing P1 pressure rows, then retain only nonzero nodes.
    /// Selected normals point OUT of the fluid; facet side is explicit. A patch
    /// must be planar with one consistent normal. Its area must later match
    /// the gas face to roundoff, so mesh mismatch is never hidden by rescaling.
    pub fn new(
        mesh: &ShellMesh,
        selections: Vec<Vec<OrientedShellFacet>>,
        budget: SurfaceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        let terms = mesh
            .nodes
            .len()
            .checked_mul(6)
            .and_then(|width| width.checked_mul(selections.len()))
            .and_then(|rows| {
                mesh.tris
                    .len()
                    .checked_mul(12)
                    .and_then(|faces| rows.checked_add(faces))
            });
        if selections.is_empty()
            || selections.len() > budget.maximum_patches
            || terms.is_none_or(|n| n > budget.maximum_preparation_terms)
        {
            return Err(Error::Invalid(
                "nonempty pressure patches exceed the preparation allowance",
            ));
        }
        let mut seen = BTreeSet::new();
        let mut patches = Vec::with_capacity(selections.len());
        for selected in selections {
            poll(&mut cancelled)?;
            if selected.is_empty() {
                return Err(Error::Invalid("empty fluid-wall patch"));
            }
            let first = mesh.facet(selected[0].triangle)?;
            let sign = orientation(selected[0].side);
            let normal = first.frame[2].map(|x| sign * x);
            let base = mesh.nodes[mesh.tris[selected[0].triangle][0]];
            for selection in &selected {
                poll(&mut cancelled)?;
                if !seen.insert(selection.triangle) {
                    return Err(Error::Invalid(
                        "an original structural facet belongs to two fluid faces",
                    ));
                }
                let facet = mesh.facet(selection.triangle)?;
                let n = facet.frame[2].map(|x| orientation(selection.side) * x);
                if !same_normal(normal, n) {
                    return Err(Error::Invalid(
                        "fluid face requires one oriented planar patch",
                    ));
                }
                for &node in &mesh.tris[selection.triangle] {
                    let delta = core::array::from_fn(|i| mesh.nodes[node][i] - base[i]);
                    let scale = delta.iter().fold(0.0_f64, |m, value| m.max(value.abs()));
                    if dot(delta, normal).abs() > ROUNDING * scale {
                        return Err(Error::Invalid("fluid-wall patch is not planar"));
                    }
                }
            }
            let port = ShellPressurePort::new(mesh, selected, budget.geometry, &mut cancelled)?;
            let mut row = Vec::new();
            for (node, weights) in port.physical_row_m2().chunks_exact(6).enumerate() {
                if node % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                let area = [weights[0], weights[1], weights[2]];
                if area != [0.0; 3] {
                    #[cfg(feature = "experimental-compressible")]
                    let scalar_area = dot(area, normal);
                    #[cfg(feature = "experimental-compressible")]
                    if !positive(scalar_area) {
                        return Err(Error::Invalid(
                            "P1 vector-traction area weight is unrepresentable",
                        ));
                    }
                    row.push(PatchWeight {
                        node,
                        normal_area_m2: area,
                        #[cfg(feature = "experimental-compressible")]
                        area_m2: scalar_area,
                    });
                }
            }
            patches.push(Patch {
                area_m2: port.area_m2(),
                normal,
                row,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            positions: mesh.nodes.clone(),
            patches,
        })
    }

    /// Original full structural node count, including unloaded nodes.
    pub fn node_count(&self) -> usize {
        self.positions.len()
    }
    /// Number of separately matched fluid-wall pressure patches.
    pub fn patch_count(&self) -> usize {
        self.patches.len()
    }

    /// Project one simultaneous accepted impulse per original fluid face.
    /// The exterior pressure acts on the opposite side. The resultant and
    /// torque are independently summed from the returned original-node loads.
    /// Supplied `origin_m` changes only the reported moment reference.
    pub fn project(
        &self,
        impulses: &[WallPressureImpulse],
        origin_m: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ShockLoad, Error> {
        self.project_impulses(
            impulses.iter().map(|value| value.inner),
            origin_m,
            cancelled,
        )
    }

    /// Project complete signed molecular wall tractions through the same
    /// original P1 rows. Tensile normal stress is retained; no positive-pressure
    /// clamp or second viscous force is applied. No interface is moved.
    #[cfg(feature = "experimental-compressible")]
    pub fn project_tractions(
        &self,
        impulses: &[WallTractionImpulse],
        origin_m: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ShockLoad, Error> {
        self.project_impulses(
            impulses.iter().map(|value| value.inner),
            origin_m,
            cancelled,
        )
    }

    /// Project complete pressure and viscous shear impulses onto the original
    /// translational P1 degrees of freedom. Scalar shape-function area weights
    /// are prepared once from the original planar pressure row, then applied
    /// to every traction component. Exterior pressure acts only along the
    /// matched normal; the tangential gas impulse is retained without clipping.
    ///
    /// The resulting nodal impulse, its moment and spatial virtual work share
    /// the same sparse projection. This remains prescribed stationary-flow
    /// loading: structural motion does not feed back into the gas wall.
    #[cfg(feature = "experimental-compressible")]
    pub fn project_vector_tractions(
        &self,
        impulses: &[WallVectorTractionImpulse],
        origin_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ShockLoad, Error> {
        poll(&mut cancelled)?;
        if impulses.len() != self.patches.len() || origin_m.iter().any(|x| !x.is_finite()) {
            return Err(Error::Invalid(
                "one wall impulse per patch and finite moment origin required",
            ));
        }
        let first = impulses
            .first()
            .ok_or(Error::Invalid("empty wall impulse list"))?;
        // Keep overlapping face contributions until all impulses at an
        // original node have been assembled. Opposing large wall loads must
        // not erase a third, smaller resolved shear before cancellation.
        let mut nodal_impulses = vec![[Dd::ZERO; 3]; self.positions.len()];
        let mut gas_impulse = [Dd::ZERO; 3];
        let mut exterior_impulse = [Dd::ZERO; 3];
        for (patch, impulse) in self.patches.iter().zip(impulses) {
            poll(&mut cancelled)?;
            if impulse.duration_s() != first.duration_s()
                || impulse.end_time_s() != first.end_time_s()
                || !same_normal(patch.normal, impulse.outward_normal())
                || (patch.area_m2 - impulse.area_m2()).abs()
                    > ROUNDING * patch.area_m2.max(impulse.area_m2())
            {
                return Err(Error::Invalid(
                    "wall/surface area, outward normal or accepted interval mismatch",
                ));
            }
            let gas = impulse.gas_traction_impulse_pa_s();
            let net = impulse.net_traction_impulse_pa_s();
            for axis in 0..3 {
                gas_impulse[axis] = gas_impulse[axis] + Dd::from_f64(gas[axis] * impulse.area_m2());
                exterior_impulse[axis] = exterior_impulse[axis]
                    - Dd::from_f64(
                        impulse.outward_normal()[axis]
                            * impulse.exterior_pressure_impulse_pa_s()
                            * impulse.area_m2(),
                    );
            }
            for (i, weight) in patch.row.iter().enumerate() {
                if i % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                for axis in 0..3 {
                    nodal_impulses[weight.node][axis] = nodal_impulses[weight.node][axis]
                        + Dd::from_f64(weight.area_m2 * net[axis]);
                }
            }
        }
        self.finish_load(
            ShockLoad {
                duration_s: first.duration_s(),
                end_time_s: first.end_time_s(),
                nodal_impulses_n_s: nodal_impulses
                    .into_iter()
                    .map(|value| value.map(Dd::to_f64))
                    .collect(),
                mean_nodal_forces_n: Vec::new(),
                gas_impulse_n_s: gas_impulse.map(Dd::to_f64),
                exterior_impulse_n_s: exterior_impulse.map(Dd::to_f64),
                resultant_impulse_n_s: [0.0; 3],
                angular_impulse_n_m_s: [0.0; 3],
            },
            origin_m,
            cancelled,
        )
    }

    fn project_impulses(
        &self,
        mut impulses: impl ExactSizeIterator<Item = WallImpulse>,
        origin_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ShockLoad, Error> {
        poll(&mut cancelled)?;
        if impulses.len() != self.patches.len() || origin_m.iter().any(|x| !x.is_finite()) {
            return Err(Error::Invalid(
                "one wall impulse per patch and finite moment origin required",
            ));
        }
        let first = impulses
            .next()
            .ok_or(Error::Invalid("empty wall impulse list"))?;
        let mut nodal_impulses = vec![[0.0; 3]; self.positions.len()];
        let mut gas_impulse = [0.0; 3];
        let mut exterior_impulse = [0.0; 3];
        for (patch, impulse) in self
            .patches
            .iter()
            .zip(std::iter::once(first).chain(impulses))
        {
            poll(&mut cancelled)?;
            if impulse.duration_s != first.duration_s
                || impulse.end_time_s != first.end_time_s
                || !same_normal(patch.normal, impulse.normal)
                || (patch.area_m2 - impulse.area_m2).abs()
                    > ROUNDING * patch.area_m2.max(impulse.area_m2)
            {
                return Err(Error::Invalid(
                    "wall/surface area, outward normal or accepted interval mismatch",
                ));
            }
            for axis in 0..3 {
                gas_impulse[axis] += impulse.normal[axis] * impulse.gas_pa_s * impulse.area_m2;
                exterior_impulse[axis] -=
                    impulse.normal[axis] * impulse.exterior_pa_s * impulse.area_m2;
            }
            for (i, weight) in patch.row.iter().enumerate() {
                if i % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                for axis in 0..3 {
                    nodal_impulses[weight.node][axis] +=
                        weight.normal_area_m2[axis] * impulse.net_pa_s();
                }
            }
        }
        self.finish_load(
            ShockLoad {
                duration_s: first.duration_s,
                end_time_s: first.end_time_s,
                nodal_impulses_n_s: nodal_impulses,
                mean_nodal_forces_n: Vec::new(),
                gas_impulse_n_s: gas_impulse,
                exterior_impulse_n_s: exterior_impulse,
                resultant_impulse_n_s: [0.0; 3],
                angular_impulse_n_m_s: [0.0; 3],
            },
            origin_m,
            cancelled,
        )
    }

    fn finish_load(
        &self,
        mut load: ShockLoad,
        origin_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ShockLoad, Error> {
        let mut forces = Vec::with_capacity(load.nodal_impulses_n_s.len());
        let mut resultant = [Dd::ZERO; 3];
        let mut moment = [Dd::ZERO; 3];
        for (i, (position, impulse)) in self
            .positions
            .iter()
            .zip(&load.nodal_impulses_n_s)
            .enumerate()
        {
            if i % 256 == 0 {
                poll(&mut cancelled)?;
            }
            let lever = core::array::from_fn(|axis| position[axis] - origin_m[axis]);
            let torque = cross(lever, *impulse);
            for axis in 0..3 {
                resultant[axis] = resultant[axis] + Dd::from_f64(impulse[axis]);
                moment[axis] = moment[axis] + Dd::from_f64(torque[axis]);
            }
            forces.push(impulse.map(|x| x / load.duration_s));
        }
        let resultant = resultant.map(Dd::to_f64);
        let moment = moment.map(Dd::to_f64);
        if forces
            .iter()
            .flatten()
            .chain(load.nodal_impulses_n_s.iter().flatten())
            .chain(&resultant)
            .chain(&moment)
            .chain(&load.gas_impulse_n_s)
            .chain(&load.exterior_impulse_n_s)
            .any(|x| !x.is_finite())
        {
            return Err(Error::Invalid(
                "structural impulse, load or moment overflow",
            ));
        }
        poll(&mut cancelled)?;
        load.mean_nodal_forces_n = forces;
        load.resultant_impulse_n_s = resultant;
        load.angular_impulse_n_m_s = moment;
        Ok(load)
    }
}

/// Original structural node order, including supported nodes. A constraint
/// owner may retain their reactions; this adapter never redistributes them.
#[derive(Clone, Debug, PartialEq)]
pub struct ShockLoad {
    /// Common accepted physical duration [s].
    pub duration_s: f64,
    /// Common accepted gas endpoint clock [s].
    pub end_time_s: f64,
    /// Integrated force on each original structural node [N s].
    pub nodal_impulses_n_s: Vec<[f64; 3]>,
    /// Corresponding interval-mean nodal forces [N], with no constraint projection.
    pub mean_nodal_forces_n: Vec<[f64; 3]>,
    /// Gas-side support impulse before exterior pressure subtraction [N s].
    pub gas_impulse_n_s: [f64; 3],
    /// Opposite impulse from the explicitly held exterior pressure [N s].
    pub exterior_impulse_n_s: [f64; 3],
    /// Sum of the returned structural nodal impulses [N s].
    pub resultant_impulse_n_s: [f64; 3],
    /// Moment of returned nodal impulses about the requested origin [N m s].
    pub angular_impulse_n_m_s: [f64; 3],
}
impl ShockLoad {
    /// Spatial virtual work of the interval-mean traction on original nodal
    /// displacements. This is not a claim about unresolved pressure-time work.
    pub fn virtual_work_j(&self, displacement_m: &[[f64; 3]]) -> Result<f64, Error> {
        if displacement_m.len() != self.mean_nodal_forces_n.len()
            || displacement_m.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(Error::Invalid(
                "finite displacement at each original structural node required",
            ));
        }
        let mut work = Dd::ZERO;
        for (force, displacement) in self.mean_nodal_forces_n.iter().zip(displacement_m) {
            for axis in 0..3 {
                work = work + Dd::from_f64(force[axis] * displacement[axis]);
            }
        }
        let work = work.to_f64();
        if !work.is_finite() {
            return Err(Error::Invalid("virtual work overflow"));
        }
        Ok(work)
    }
}

const ROUNDING: f64 = 1024.0 * f64::EPSILON;
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn unit(x: [f64; 3]) -> bool {
    x.iter().all(|x| x.is_finite()) && (dot(x, x) - 1.0).abs() <= ROUNDING
}
fn same_normal(a: [f64; 3], b: [f64; 3]) -> bool {
    a.into_iter().zip(b).all(|(a, b)| (a - b).abs() <= ROUNDING)
}
fn orientation(side: ShellFace) -> f64 {
    match side {
        ShellFace::Positive => 1.0,
        ShellFace::Negative => -1.0,
    }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
