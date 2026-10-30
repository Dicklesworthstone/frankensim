//! Three-dimensional traces using the original mixture HLLC/HLL owner.

use crate::multicomponent::{
    Cell, Error, FaceFlux, MAX_COMPONENTS, MAX_SPATIAL_COMPONENTS, MixtureGas, NormalCell,
    Primitive, poll, sum,
};

/// Borrowed conservative mixture trace in an orthonormal three-dimensional
/// face frame. No species is treated as an independently advected passive scalar.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialCell<'a> {
    /// Original canonical species partial mass densities [kg/m³].
    pub partial_densities: &'a [f64],
    /// Momentum density along the oriented unit face normal [kg/(m² s)].
    pub normal_momentum: f64,
    /// Momentum densities along two orthonormal tangents [kg/(m² s)].
    pub transverse_momentum: [f64; 2],
    /// Full total energy less the original species references [J/m³].
    pub energy: f64,
}
impl SpatialCell<'_> {
    /// Mixture density is the compensated sum of its species partial densities.
    pub fn density(self) -> f64 {
        sum(self.partial_densities.iter().copied())
    }
}
impl<'a> From<Cell<'a>> for SpatialCell<'a> {
    fn from(cell: Cell<'a>) -> Self {
        Self {
            partial_densities: cell.partial_densities,
            normal_momentum: cell.momentum,
            transverse_momentum: [0.0; 2],
            energy: cell.energy,
        }
    }
}
impl<'a> From<NormalCell<'a>> for SpatialCell<'a> {
    fn from(cell: NormalCell<'a>) -> Self {
        Self {
            partial_densities: cell.partial_densities,
            normal_momentum: cell.normal_momentum,
            transverse_momentum: [cell.transverse_momentum, 0.0],
            energy: cell.energy,
        }
    }
}

/// Reconstructed mixture trace, retaining original species order and calorics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialPrimitive<'a> {
    /// Current species partial densities [kg/m³], including exact-zero species.
    pub partial_densities: &'a [f64],
    /// Normal velocity [m/s].
    pub normal_velocity: f64,
    /// Velocities along both face tangents [m/s].
    pub transverse_velocity: [f64; 2],
    /// Absolute pressure [Pa], evaluated through the original mixture source.
    pub pressure: f64,
}

impl MixtureGas {
    /// Construct a three-dimensional conservative trace and its original source
    /// thermodynamics without allocating. Output has exactly `species_count()+4`
    /// entries: species, normal momentum, two transverse momenta, full referenced
    /// total energy. The returned primitive velocity is normal.
    ///
    /// All kinetic components are included before checking the declared caloric
    /// precision. Output is scratch until success; unresolved sensible energy at
    /// extreme velocity is refused, never reconstructed by clipping temperature.
    pub fn spatial_conserved_into(
        &self,
        input: SpatialPrimitive<'_>,
        out: &mut [f64],
    ) -> Result<Primitive, Error> {
        let n = self.species_count();
        if out.len() != n + 4 || !input.transverse_velocity.iter().all(|v| v.is_finite()) {
            return Err(Error::Invalid(
                "spatial mixture output axis or transverse velocity",
            ));
        }
        let mut normal = [0.0; MAX_COMPONENTS];
        let primitive = self.state_at(
            input.partial_densities,
            input.normal_velocity,
            input.pressure,
            &mut normal[..n + 2],
        )?;
        out[..n + 1].copy_from_slice(&normal[..n + 1]);
        out[n + 3] = normal[n + 1];
        for (axis, velocity) in input.transverse_velocity.into_iter().enumerate() {
            out[n + 1 + axis] = primitive.density * velocity;
            out[n + 3] += primitive.density * (0.5 * velocity * velocity);
        }
        let reduced = self.spatial_reduction(from_slice(out, n))?;
        let recovered =
            reduced.energy - 0.5 * reduced.momentum * (reduced.momentum / primitive.density);
        let (thermal, _) = self.caloric(input.partial_densities, primitive.temperature)?;
        let tolerance = primitive.density * self.domain().absolute_energy_tolerance
            + self.domain().relative_energy_tolerance * thermal.abs();
        if !self.admissible_cell(reduced) || (recovered - thermal).abs() > tolerance {
            return Err(Error::Invalid(
                "spatial mixture loses declared caloric precision",
            ));
        }
        Ok(primitive)
    }

    /// Decode all three kinetic components through the original evolving-mixture
    /// caloric inversion. The returned primitive velocity is normal.
    pub fn spatial_primitive(
        &self,
        input: SpatialCell<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        poll(&mut cancelled)?;
        self.primitive(self.spatial_reduction(input)?, &mut cancelled)
    }

    /// Shared species/three-momentum/energy Riemann flux with full kinetic-energy
    /// admission for both HLLC stars and the bounded HLL fallback. Output follows
    /// the constructor's `species_count()+4` layout and is scratch until success.
    pub fn hllc_spatial_flux(
        &self,
        left: SpatialCell<'_>,
        right: SpatialCell<'_>,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let n = self.species_count();
        if out.len() != n + 4 {
            return Err(Error::Invalid("spatial mixture flux output axis"));
        }
        let pl = self.spatial_primitive(left, &mut cancelled)?;
        let pr = self.spatial_primitive(right, &mut cancelled)?;
        let mut l = [0.0; MAX_SPATIAL_COMPONENTS];
        let mut r = [0.0; MAX_SPATIAL_COMPONENTS];
        write_cell(left, &mut l[..n + 4]);
        write_cell(right, &mut r[..n + 4]);
        self.interface_components::<2>(&l[..n + 4], pl, &r[..n + 4], pr, out, &mut cancelled)
    }

    /// Reuse original source thermodynamics computed while admitting reconstructed
    /// traces, avoiding a second temperature inversion at each spatial face.
    pub fn hllc_spatial_flux_from_primitive(
        &self,
        left: SpatialPrimitive<'_>,
        right: SpatialPrimitive<'_>,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(&mut cancelled)?;
        let n = self.species_count();
        if out.len() != n + 4 {
            return Err(Error::Invalid("spatial mixture flux output axis"));
        }
        let mut l = [0.0; MAX_SPATIAL_COMPONENTS];
        let mut r = [0.0; MAX_SPATIAL_COMPONENTS];
        let pl = self.spatial_conserved_into(left, &mut l[..n + 4])?;
        poll(&mut cancelled)?;
        let pr = self.spatial_conserved_into(right, &mut r[..n + 4])?;
        self.interface_components::<2>(&l[..n + 4], pl, &r[..n + 4], pr, out, &mut cancelled)
    }

    /// Original composition-local shock/rarefaction pressure at a stationary
    /// adiabatic slip wall. The outward normal is -1 or +1. All species, both
    /// transverse momenta and total energy have exactly zero wall flux.
    pub fn stationary_wall_spatial_flux(
        &self,
        inside: SpatialCell<'_>,
        outward_normal: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let p = self.spatial_primitive(inside, &mut cancelled)?;
        self.spatial_wall(inside, p, outward_normal, out, &mut cancelled)
    }

    /// Source-admitted reconstructed trace through the same physical wall owner.
    pub fn stationary_wall_spatial_flux_from_primitive(
        &self,
        inside: SpatialPrimitive<'_>,
        outward_normal: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(&mut cancelled)?;
        let n = self.species_count();
        let mut q = [0.0; MAX_SPATIAL_COMPONENTS];
        let p = self.spatial_conserved_into(inside, &mut q[..n + 4])?;
        self.spatial_wall(
            from_slice(&q[..n + 4], n),
            p,
            outward_normal,
            out,
            &mut cancelled,
        )
    }

    fn spatial_wall(
        &self,
        inside: SpatialCell<'_>,
        primitive: Primitive,
        outward_normal: f64,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let n = self.species_count();
        if out.len() != n + 4 {
            return Err(Error::Invalid("spatial mixture wall output axis"));
        }
        let mut normal = [0.0; MAX_COMPONENTS];
        let report = self.wall_flux(
            self.spatial_reduction(inside)?,
            primitive,
            outward_normal,
            &mut normal[..n + 2],
            cancelled,
        )?;
        out.fill(0.0);
        out[n] = normal[n];
        Ok(report)
    }

    pub(in crate::multicomponent) fn admissible_spatial_cell(&self, cell: SpatialCell<'_>) -> bool {
        self.spatial_reduction(cell)
            .is_ok_and(|q| self.admissible_cell(q))
    }

    pub(in crate::multicomponent) fn spatial_reduction<'a>(
        &self,
        cell: SpatialCell<'a>,
    ) -> Result<Cell<'a>, Error> {
        let (rho, _) = self.admit_densities(cell.partial_densities)?;
        if ![cell.normal_momentum, cell.energy]
            .iter()
            .chain(&cell.transverse_momentum)
            .all(|x| x.is_finite())
        {
            return Err(Error::Invalid("non-finite spatial mixture trace"));
        }
        let mut energy = cell.energy;
        for momentum in cell.transverse_momentum {
            energy -= 0.5 * momentum * (momentum / rho);
        }
        if !energy.is_finite() {
            return Err(Error::Invalid(
                "unrepresentable spatial mixture kinetic energy",
            ));
        }
        Ok(Cell {
            partial_densities: cell.partial_densities,
            momentum: cell.normal_momentum,
            energy,
        })
    }
}

fn from_slice(q: &[f64], n: usize) -> SpatialCell<'_> {
    SpatialCell {
        partial_densities: &q[..n],
        normal_momentum: q[n],
        transverse_momentum: [q[n + 1], q[n + 2]],
        energy: q[n + 3],
    }
}
fn write_cell(cell: SpatialCell<'_>, out: &mut [f64]) {
    let n = cell.partial_densities.len();
    out[..n].copy_from_slice(cell.partial_densities);
    out[n] = cell.normal_momentum;
    out[n + 1..n + 3].copy_from_slice(&cell.transverse_momentum);
    out[n + 3] = cell.energy;
}

#[cfg(test)]
mod tests;
