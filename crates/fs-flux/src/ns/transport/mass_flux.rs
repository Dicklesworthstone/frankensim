//! Compatible mass flux for conservative P0 density and BDM1 momentum.
//!
//! The directional donor density multiplies the SAME original BDM face flow
//! used by density transport. A reversing face has a piecewise linear mass
//! trace. Its moments through degree two are essential: BDM1 kinetic energy
//! is quadratic, so keeping only the mean or a BDM1 trace loses the discrete
//! integration-by-parts identity needed by the momentum energy balance.
//! With bounded linear density, each directional donor is the exact volume-
//! flux-weighted mean of its limited affine density. It still multiplies the
//! original affine normal velocity in this momentum owner. Thus the accepted
//! directional mass and kinetic moments remain compatible, without claiming
//! a pointwise quadratic mass trace or higher spatial order for momentum.
//!
//! On each reference triangle, a vector P2 polynomial matches all nine face
//! moments and has constant divergence. Its one remaining divergence-free,
//! zero-normal bubble minimizes physical L2 distance to rho_cell * velocity.
//! This is an explicit local BDM2 reconstruction, without a global solve or
//! per-cell matrix factorization. The geometry uses the contravariant Piola
//! map; constant divergence is the original P0 mass balance divided by area.

use super::density::{DensityTransportCell, DensityTransportError, PreparedDensityTransport};
use super::dimensional_rate;
use crate::bdm::tri_quad;
use crate::ns::convection_edge_quadrature;
use fs_math::dd::Dd;

/// Native BDM1 momentum transported by this same conservative mass flux.
pub mod momentum;

type Result<T> = core::result::Result<T, DensityTransportError>;

/// One original sign-split face site. Both rates ALREADY include its line
/// quadrature weight. Positive means along the original global edge normal.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DensityMassFluxPoint {
    pub point_m: [f64; 2],
    /// Signed volume rate per depth [m²/s], independent of donor density.
    pub volume_rate_m2_s: f64,
    /// Signed mass rate per depth [kg/(m s)].
    pub mass_rate_kg_m_s: f64,
    /// Affine interpolation fraction between the original two boundary sites.
    pub boundary_fraction: f64,
}

struct Face {
    points: [DensityMassFluxPoint; 4],
    count: usize,
    directional: [f64; 2],
    // Integrals against shifted Legendre polynomials 1, 2s-1, 6s²-6s+1.
    moments: [f64; 3],
}

struct Cell {
    origin: [f64; 2],
    length: f64,
    shape: [[f64; 2]; 2],
    determinant: f64,
    flux_scale: f64,
    // Reference vector P2 in [1, x, y, x², xy, y²], divided by flux_scale.
    coefficients: [[f64; 6]; 2],
    divergence: f64,
}

/// Frozen original face streams and their conservative cell mass flux.
/// Storage and preparation are linear in cells+edges. Sampling allocates
/// nothing. This reconstructs transport; it neither advances velocity nor
/// supplies a pressure, constitutive law, EOS, or new continuum-order claim.
pub struct PreparedDensityMassFlux<'m> {
    mesh: &'m crate::TriMesh,
    faces: Vec<Face>,
    cells: Vec<Cell>,
}

impl<'m> PreparedDensityTransport<'m> {
    /// Reconstruct the mass flux of these exact frozen volume faces and donor
    /// states. Use the SAME states in `advance_into` and the momentum step.
    /// Input quantities per unit mass are validated but do not affect flux.
    /// No timestep or physical depth multiplies the returned per-depth rates.
    ///
    /// Face moments through degree two reproduce the piecewise donor flux;
    /// the constant cell divergence equals net outward mass rate / cell area.
    /// Consequently the Euler density increment and divergence agree up to
    /// the represented-state roundoff already reported by density transport.
    /// Both signs of a reversing face are retained, including zero net volume.
    ///
    /// `max_storage_bytes` bounds the retained vectors, excluding this original
    /// plan and caller inputs. Bounded linear density also admits its temporary
    /// vertex/directional reconstruction arrays within this same allowance.
    /// Errors and per-face/per-cell cancellation publish no partial operator.
    pub fn prepare_mass_flux<const N: usize>(
        &self,
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        max_storage_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedDensityMassFlux<'m>> {
        poll(&mut cancelled)?;
        let required = mass_flux_storage(self.mesh)?;
        let donors = self.reconstructed_mass_flux_donors(
            cells,
            inflow,
            required,
            max_storage_bytes,
            &mut cancelled,
        )?;
        self.prepare_mass_flux_with_donors(
            cells,
            inflow,
            donors.as_deref(),
            max_storage_bytes,
            &mut cancelled,
        )
    }

    // Time integration passes the actual mean of its limited directional
    // donors. Reconstructing an averaged cell state would change the mass flux
    // because the spatial limiter is nonlinear.
    pub(super) fn prepare_mass_flux_with_donors<const N: usize>(
        &self,
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        directional_donors: Option<&[[f64; 2]]>,
        max_storage_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedDensityMassFlux<'m>> {
        poll(&mut cancelled)?;
        let mesh = self.mesh;
        if cells.len() != mesh.tris.len()
            || inflow.len() != mesh.edges.len()
            || directional_donors.is_some_and(|donors| donors.len() != mesh.edges.len())
        {
            return Err(invalid(
                "original cell and boundary-state dimensions required",
            ));
        }
        let required = mass_flux_storage(mesh)?;
        if required > max_storage_bytes {
            return Err(DensityTransportError::Storage {
                required,
                maximum: max_storage_bytes,
            });
        }
        for state in cells {
            poll(&mut cancelled)?;
            admit_state(state)?;
        }
        let mut faces = Vec::new();
        faces
            .try_reserve_exact(mesh.edges.len())
            .map_err(|_| DensityTransportError::Allocation)?;
        for (index, edge) in mesh.edges.iter().enumerate() {
            poll(&mut cancelled)?;
            let (first, second) = edge.tris;
            if second != usize::MAX && inflow[index].is_some() {
                return Err(invalid("inflow states belong only to boundary edges"));
            }
            if let Some(state) = &inflow[index] {
                admit_state(state)?;
            }
            let rates = self.edge_volume_rates()[index];
            let original_density = [
                cells[first].density_kg_m3,
                if second != usize::MAX {
                    cells[second].density_kg_m3
                } else if rates.against_normal_m2_s() > 0.0 {
                    inflow[index]
                        .as_ref()
                        .ok_or(invalid("every incoming stream needs a donor density"))?
                        .density_kg_m3
                } else {
                    cells[first].density_kg_m3
                },
            ];
            let density = directional_donors.map_or(original_density, |donors| donors[index]);
            if density.iter().any(|rho| !rho.is_finite() || *rho <= 0.0) {
                return Err(invalid(
                    "positive reconstructed directional donor densities required",
                ));
            }
            let directional = [
                product(density[0], rates.along_normal_m2_s())?,
                product(density[1], rates.against_normal_m2_s())?,
            ];
            let a = *mesh
                .verts
                .get(edge.verts.0)
                .ok_or(invalid("original edge vertex"))?;
            let b = *mesh
                .verts
                .get(edge.verts.1)
                .ok_or(invalid("original edge vertex"))?;
            if a.into_iter().chain(b).any(|x| !x.is_finite()) {
                return Err(invalid("finite original edge coordinates required"));
            }
            let moments = self.velocity_moments[index];
            let scale = moments[0].abs().max(moments[1].abs());
            let mut face = Face {
                points: [DensityMassFluxPoint::default(); 4],
                count: 0,
                directional,
                moments: [0.0; 3],
            };
            let mut sums = [Dd::ZERO; 3];
            for (q, weight, fraction) in convection_edge_quadrature([0.0, 0.0], [1.0, 0.0], moments)
            {
                let s = q[0];
                let normal = if scale == 0.0 {
                    0.0
                } else {
                    moments[0] / scale + 12.0 * (moments[1] / scale) * (s - 0.5)
                };
                let direction = usize::from(normal < 0.0);
                let sign = if direction == 0 { 1.0 } else { -1.0 };
                let volume = if normal == 0.0 {
                    0.0
                } else {
                    dimensional_rate(weight * normal.abs(), scale, edge.len).ok_or(
                        DensityTransportError::Arithmetic("mass-flux volume quadrature"),
                    )?
                };
                let mass = sign * product(density[direction], volume)?;
                let point = core::array::from_fn(|axis| a[axis] + s * (b[axis] - a[axis]));
                if point.iter().any(|x| !x.is_finite()) {
                    return Err(DensityTransportError::Arithmetic(
                        "mass-flux quadrature position",
                    ));
                }
                face.points[face.count] = DensityMassFluxPoint {
                    point_m: point,
                    volume_rate_m2_s: sign * volume,
                    mass_rate_kg_m_s: mass,
                    boundary_fraction: fraction,
                };
                face.count += 1;
                for (degree, polynomial) in [1.0, 2.0 * s - 1.0, (6.0 * s - 6.0) * s + 1.0]
                    .into_iter()
                    .enumerate()
                {
                    sums[degree] = sums[degree] + Dd::from_f64(mass) * Dd::from_f64(polynomial);
                }
            }
            face.moments = sums.map(|value| value.to_f64());
            // The known mean uses the original directional rate, independent
            // of a second quadrature sum. A single donor gives a truly affine
            // trace: retain its exact zero quadratic coefficient.
            if rates.along_normal_m2_s() == 0.0
                || rates.against_normal_m2_s() == 0.0
                || density[0] == density[1]
            {
                let rho = if rates.along_normal_m2_s() == 0.0 {
                    density[1]
                } else {
                    density[0]
                };
                face.moments = [
                    signed_product(rho, rates.net_m2_s())?,
                    signed_product(rho, 2.0 * edge.len * moments[1])?,
                    0.0,
                ];
            } else {
                face.moments[0] =
                    (Dd::from_f64(directional[0]) - Dd::from_f64(directional[1])).to_f64();
            }
            if face.moments.iter().any(|x| !x.is_finite()) {
                return Err(DensityTransportError::Arithmetic("mass-flux face moments"));
            }
            faces.push(face);
        }
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(cells.len())
            .map_err(|_| DensityTransportError::Allocation)?;
        for (cell, state) in cells.iter().enumerate() {
            poll(&mut cancelled)?;
            let tri = mesh.tris[cell];
            let mut p = [[0.0; 2]; 3];
            for k in 0..3 {
                p[k] = *mesh
                    .verts
                    .get(tri[k])
                    .ok_or(invalid("original triangle vertex"))?;
            }
            let map = [
                [p[1][0] - p[0][0], p[2][0] - p[0][0]],
                [p[1][1] - p[0][1], p[2][1] - p[0][1]],
            ];
            let length = map.iter().flatten().fold(0.0_f64, |a, b| a.max(b.abs()));
            if !length.is_finite() || length <= 0.0 {
                return Err(invalid("finite nondegenerate reference triangle required"));
            }
            let shape = map.map(|row| row.map(|x| x / length));
            let determinant = shape[0][0] * shape[1][1] - shape[0][1] * shape[1][0];
            let area = 0.5 * (map[0][0] * map[1][1] - map[0][1] * map[1][0]);
            if !determinant.is_finite()
                || determinant <= 0.0
                || !area.is_finite()
                || area <= 0.0
                || (area - mesh.areas[cell]).abs()
                    > 64.0 * f64::EPSILON * area.max(mesh.areas[cell])
            {
                return Err(invalid("original CCW triangle area required"));
            }
            let mut mass = [[0.0; 3]; 3];
            let mut target = [[0.0; 3]; 3];
            let mut signs = [0.0; 3];
            let mut orientation = [0.0; 3];
            let mut outward = Dd::ZERO;
            let mut flux_scale = 0.0_f64;
            for k in 0..3 {
                let (edge_id, sign) = mesh.tri_edges[cell][k];
                let edge = &mesh.edges[edge_id];
                let (first, second) = (tri[(k + 1) % 3], tri[(k + 2) % 3]);
                if edge.verts != (first.min(second), first.max(second)) {
                    return Err(invalid("original opposite-vertex edge topology required"));
                }
                let tangent = [
                    p[(k + 2) % 3][0] - p[(k + 1) % 3][0],
                    p[(k + 2) % 3][1] - p[(k + 1) % 3][1],
                ];
                let normal_line = [tangent[1], -tangent[0]];
                for axis in 0..2 {
                    if !edge.normal[axis].is_finite()
                        || (sign * edge.len * edge.normal[axis] - normal_line[axis]).abs()
                            > 64.0 * f64::EPSILON * edge.len
                    {
                        return Err(invalid(
                            "original outward unit normals and line measures required",
                        ));
                    }
                }
                signs[k] = sign;
                orientation[k] = if edge.verts.0 == first { 1.0 } else { -1.0 };
                mass[k] = faces[edge_id].moments;
                target[k] = [
                    signed_product(
                        state.density_kg_m3,
                        self.edge_volume_rates()[edge_id].net_m2_s(),
                    )?,
                    signed_product(
                        state.density_kg_m3,
                        2.0 * edge.len * self.velocity_moments[edge_id][1],
                    )?,
                    0.0,
                ];
                outward = outward + Dd::from_f64(sign * mass[k][0]);
                for value in mass[k].into_iter().chain(target[k]) {
                    flux_scale = flux_scale.max(value.abs());
                }
            }
            let mut coefficients = [[0.0; 6]; 2];
            if flux_scale > 0.0 {
                let traces = core::array::from_fn(|k| {
                    trace(mass[k].map(|x| x / flux_scale), signs[k], orientation[k])
                });
                let target_traces = core::array::from_fn(|k| {
                    trace(target[k].map(|x| x / flux_scale), signs[k], orientation[k])
                });
                coefficients = reconstruct(traces);
                let target = reconstruct(target_traces);
                // The zero-normal bubble B=(-d_y b,d_x b), b=xy(1-x-y).
                // A^T A supplies the PHYSICAL L2 metric; using the reference
                // Euclidean norm would depend on the local vertex numbering.
                let mut numerator = Dd::ZERO;
                let mut denominator = Dd::ZERO;
                for (q, weight) in tri_quad([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], 0.5) {
                    let bubble = [
                        -q[0] + q[0] * q[0] + 2.0 * q[0] * q[1],
                        q[1] - 2.0 * q[0] * q[1] - q[1] * q[1],
                    ];
                    let difference = core::array::from_fn(|axis| {
                        evaluate(&coefficients[axis], q) - evaluate(&target[axis], q)
                    });
                    let b = multiply(shape, bubble);
                    let d = multiply(shape, difference);
                    for axis in 0..2 {
                        numerator = numerator
                            + Dd::from_f64(weight) * Dd::from_f64(b[axis]) * Dd::from_f64(d[axis]);
                        denominator = denominator
                            + Dd::from_f64(weight) * Dd::from_f64(b[axis]) * Dd::from_f64(b[axis]);
                    }
                }
                let amplitude = -numerator.to_f64() / denominator.to_f64();
                coefficients[0][1] -= amplitude;
                coefficients[0][3] = amplitude;
                coefficients[1][2] += amplitude;
                coefficients[1][4] = -2.0 * amplitude;
                coefficients[1][5] -= amplitude;
                coefficients[0][4] = -2.0 * coefficients[1][5];
            }
            let divergence = outward.to_f64() / area;
            if !divergence.is_finite() || coefficients.iter().flatten().any(|x| !x.is_finite()) {
                return Err(DensityTransportError::Arithmetic(
                    "mass-flux reconstruction",
                ));
            }
            prepared.push(Cell {
                origin: p[0],
                length,
                shape,
                determinant,
                flux_scale,
                coefficients,
                divergence,
            });
        }
        poll(&mut cancelled)?;
        Ok(PreparedDensityMassFlux {
            mesh,
            faces,
            cells: prepared,
        })
    }
}

impl PreparedDensityMassFlux<'_> {
    /// Weighted original face quadrature, including separate flow directions.
    pub fn edge_quadrature(&self, edge: usize) -> Option<&[DensityMassFluxPoint]> {
        self.faces.get(edge).map(|face| &face.points[..face.count])
    }
    /// Nonnegative directional mass rates [kg/(m s)], in original normal order.
    pub fn edge_directional_mass_rates(&self, edge: usize) -> Option<[f64; 2]> {
        self.faces.get(edge).map(|face| face.directional)
    }
    /// Constant divergence [kg/(m³ s)] of the local conservative reconstruction.
    pub fn cell_divergence(&self, cell: usize) -> Option<f64> {
        self.cells.get(cell).map(|cell| cell.divergence)
    }
    /// Original mesh identity, used by native consumers to reject a foreign map.
    pub fn mesh(&self) -> &crate::TriMesh {
        self.mesh
    }

    /// Local polynomial mass flux [kg/(m² s)] at a finite physical point.
    /// `cell` explicitly selects a one-sided volume polynomial; it is not a
    /// spatial search or an average across material interfaces. Callers own
    /// support membership. Reference mapping avoids global-coordinate powers.
    pub fn cell_value(&self, cell: usize, point_m: [f64; 2]) -> Result<[f64; 2]> {
        let cell = self
            .cells
            .get(cell)
            .ok_or(invalid("mass-flux cell index"))?;
        if point_m.iter().any(|x| !x.is_finite()) {
            return Err(invalid("finite mass-flux point required"));
        }
        let delta = core::array::from_fn::<_, 2, _>(|axis| {
            (point_m[axis] - cell.origin[axis]) / cell.length
        });
        let h = cell.shape;
        let q = [
            (h[1][1] * delta[0] - h[0][1] * delta[1]) / cell.determinant,
            (h[0][0] * delta[1] - h[1][0] * delta[0]) / cell.determinant,
        ];
        let reference = cell.coefficients.map(|c| evaluate(&c, q));
        let mapped = multiply(h, reference);
        let denominator = cell.length * cell.determinant;
        Ok([
            scaled_quotient(mapped[0], cell.flux_scale, denominator)?,
            scaled_quotient(mapped[1], cell.flux_scale, denominator)?,
        ])
    }
}

fn trace(m: [f64; 3], sign: f64, orientation: f64) -> [f64; 3] {
    [
        sign * (m[0] - 3.0 * orientation * m[1] + 5.0 * m[2]),
        sign * (6.0 * orientation * m[1] - 30.0 * m[2]),
        sign * 30.0 * m[2],
    ]
}

// Reference outward line normals on edges opposite (0,0),(1,0),(0,1)
// are (1,1),(-1,0),(0,-1). The two nonconstant divergence terms vanish.
// The free bubble amplitude is initially zero and is selected in physical L2.
fn reconstruct(g: [[f64; 3]; 3]) -> [[f64; 6]; 2] {
    let a0 = -(g[1][0] + g[1][1] + g[1][2]);
    let a2 = g[1][1] + 2.0 * g[1][2];
    let a5 = -g[1][2];
    let b0 = -g[2][0];
    let b1 = -g[2][1];
    let b3 = -g[2][2];
    let r0 = g[0][0] - a0 - b0 - b1 - b3;
    let r1 = g[0][1] - a2 + b1 + 2.0 * b3;
    let r2 = (g[0][2] - a5 - b3) / 3.0;
    [
        [a0, r0, a2, 0.0, -2.0 * r2, a5],
        [b0, b1, r0 + r1 + 2.0 * r2, b3, 0.0, r2],
    ]
}

fn evaluate(c: &[f64; 6], q: [f64; 2]) -> f64 {
    c[0] + q[0] * (c[1] + q[0] * c[3] + q[1] * c[4]) + q[1] * (c[2] + q[1] * c[5])
}
fn multiply(a: [[f64; 2]; 2], x: [f64; 2]) -> [f64; 2] {
    [
        a[0][0] * x[0] + a[0][1] * x[1],
        a[1][0] * x[0] + a[1][1] * x[1],
    ]
}
fn scaled_quotient(value: f64, scale: f64, divisor: f64) -> Result<f64> {
    if !value.is_finite() || !divisor.is_finite() || divisor <= 0.0 {
        return Err(DensityTransportError::Arithmetic("mass-flux Piola map"));
    }
    if value == 0.0 || scale == 0.0 {
        return Ok(0.0);
    }
    for result in [
        (value * scale) / divisor,
        (value / divisor) * scale,
        value * (scale / divisor),
    ] {
        if result.is_finite() && result != 0.0 {
            return Ok(result);
        }
    }
    Err(DensityTransportError::Arithmetic("mass-flux Piola range"))
}
fn product(a: f64, b: f64) -> Result<f64> {
    if b == 0.0 {
        return Ok(0.0);
    }
    let value = a * b;
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(DensityTransportError::Arithmetic("mass-flux rate"))
    }
}
fn mass_flux_storage(mesh: &crate::TriMesh) -> Result<usize> {
    mesh.tris
        .len()
        .checked_mul(size_of::<Cell>())
        .and_then(|cells| {
            mesh.edges
                .len()
                .checked_mul(size_of::<Face>())
                .and_then(|faces| cells.checked_add(faces))
        })
        .ok_or(invalid("mass-flux storage size overflow"))
}
fn signed_product(a: f64, b: f64) -> Result<f64> {
    product(a, b.abs()).map(|v| v.copysign(b))
}
fn admit_state<const N: usize>(state: &DensityTransportCell<N>) -> Result<()> {
    if !state.density_kg_m3.is_finite()
        || state.density_kg_m3 <= 0.0
        || state.specific.iter().any(|x| !x.is_finite())
    {
        Err(invalid(
            "positive finite density and finite carried state required",
        ))
    } else {
        Ok(())
    }
}
fn invalid(message: &'static str) -> DensityTransportError {
    DensityTransportError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(DensityTransportError::Cancelled)
    } else {
        Ok(())
    }
}
