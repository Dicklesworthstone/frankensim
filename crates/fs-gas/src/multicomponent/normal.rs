//! Planar traces of the original evolving-species caloric and Riemann owners.
//!
//! Spatial consumers rotate the two momentum components into an orthonormal
//! face frame. Every species, both momenta and full referenced total energy
//! then use one shared HLLC/HLL solve. In particular, transverse shear mixing
//! contributes to the HLL state's heat before the source-domain admission.

use super::{
    Cell, Error, FaceFlux, MAX_COMPONENTS, MAX_NORMAL_COMPONENTS, MixtureGas, Primitive, poll,
};

/// Borrowed conservative trace in a face's normal/tangent frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalCell<'a> {
    /// Partial species mass densities [kg/m³], in the original canonical order.
    pub partial_densities: &'a [f64],
    /// Momentum density along the oriented unit normal [kg/(m² s)].
    pub normal_momentum: f64,
    /// Momentum density along its unit tangent [kg/(m² s)].
    pub transverse_momentum: f64,
    /// Total energy including BOTH kinetic terms, less species references [J/m³].
    pub energy: f64,
}

/// Reconstructed trace; callers never provide a cached caloric model or sound speed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalPrimitive<'a> {
    /// Current partial densities [kg/m³], including exact-zero species.
    pub partial_densities: &'a [f64],
    /// Normal velocity [m/s].
    pub normal_velocity: f64,
    /// Transverse velocity [m/s].
    pub transverse_velocity: f64,
    /// Absolute pressure [Pa].
    pub pressure: f64,
}

impl<'a> From<Cell<'a>> for NormalCell<'a> {
    fn from(cell: Cell<'a>) -> Self {
        Self {
            partial_densities: cell.partial_densities,
            normal_momentum: cell.momentum,
            transverse_momentum: 0.0,
            energy: cell.energy,
        }
    }
}

impl MixtureGas {
    /// Construct an allocation-free planar conservative trace and its source
    /// thermodynamics. Output layout is species, normal momentum, transverse
    /// momentum, total referenced energy, with exactly `species_count()+3` entries.
    /// The returned primitive velocity is normal. Output is usable only on success.
    pub fn normal_conserved_into(
        &self,
        input: NormalPrimitive<'_>,
        out: &mut [f64],
    ) -> Result<Primitive, Error> {
        let n = self.species_count();
        if out.len() != n + 3 || !input.transverse_velocity.is_finite() {
            return Err(Error::Invalid(
                "planar mixture output axis or transverse velocity",
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
        out[n + 1] = primitive.density * input.transverse_velocity;
        out[n + 2] = normal[n + 1]
            + primitive.density * (0.5 * input.transverse_velocity * input.transverse_velocity);
        let cell = from_slice(out, n);
        let reduced = self.normal_reduction(cell)?;
        let recovered =
            reduced.energy - 0.5 * reduced.momentum * (reduced.momentum / primitive.density);
        let (thermal, _) = self.caloric(input.partial_densities, primitive.temperature)?;
        let tolerance = primitive.density * self.domain().absolute_energy_tolerance
            + self.domain().relative_energy_tolerance * thermal.abs();
        if !self.admissible_normal_cell(cell) || (recovered - thermal).abs() > tolerance {
            return Err(Error::Invalid(
                "planar mixture loses declared caloric precision",
            ));
        }
        Ok(primitive)
    }

    /// Decode through the existing variable-composition caloric inversion,
    /// subtracting both kinetic components. Returned velocity is normal.
    pub fn normal_primitive(
        &self,
        input: NormalCell<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        poll(&mut cancelled)?;
        self.primitive(self.normal_reduction(input)?, &mut cancelled)
    }

    /// One shared conservative species/normal/transverse/energy flux. Successful
    /// evaluation allocates nothing. `out` is scratch until success and has
    /// exactly `species_count()+3` entries in the constructor's layout.
    pub fn hllc_normal_flux(
        &self,
        left: NormalCell<'_>,
        right: NormalCell<'_>,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let n = self.species_count();
        if out.len() != n + 3 {
            return Err(Error::Invalid("planar mixture flux output axis"));
        }
        let pl = self.normal_primitive(left, &mut cancelled)?;
        let pr = self.normal_primitive(right, &mut cancelled)?;
        let mut l = [0.0; MAX_NORMAL_COMPONENTS];
        let mut r = [0.0; MAX_NORMAL_COMPONENTS];
        write_cell(left, &mut l[..n + 3]);
        write_cell(right, &mut r[..n + 3]);
        self.interface_components(&l[..n + 3], pl, &r[..n + 3], pr, out, true, &mut cancelled)
    }

    /// Source-admit reconstructed traces and reuse their computed thermodynamics
    /// directly in the original flux kernel, avoiding a second energy inversion.
    pub fn hllc_normal_flux_from_primitive(
        &self,
        left: NormalPrimitive<'_>,
        right: NormalPrimitive<'_>,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(&mut cancelled)?;
        let n = self.species_count();
        if out.len() != n + 3 {
            return Err(Error::Invalid("planar mixture flux output axis"));
        }
        let mut l = [0.0; MAX_NORMAL_COMPONENTS];
        let mut r = [0.0; MAX_NORMAL_COMPONENTS];
        let pl = self.normal_conserved_into(left, &mut l[..n + 3])?;
        poll(&mut cancelled)?;
        let pr = self.normal_conserved_into(right, &mut r[..n + 3])?;
        self.interface_components(&l[..n + 3], pl, &r[..n + 3], pr, out, true, &mut cancelled)
    }

    /// Original caloric shock/rarefaction wall pressure at a stationary slip
    /// wall. Outward normal is -1 or +1; species, tangential momentum and energy
    /// fluxes are exactly zero. Source composition is the incident cell's.
    pub fn stationary_wall_normal_flux(
        &self,
        inside: NormalCell<'_>,
        outward_normal: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let p = self.normal_primitive(inside, &mut cancelled)?;
        self.normal_wall(inside, p, outward_normal, out, &mut cancelled)
    }

    /// Source-admitted reconstructed trace through that same wall-pressure owner.
    pub fn stationary_wall_normal_flux_from_primitive(
        &self,
        inside: NormalPrimitive<'_>,
        outward_normal: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(&mut cancelled)?;
        let n = self.species_count();
        let mut q = [0.0; MAX_NORMAL_COMPONENTS];
        let p = self.normal_conserved_into(inside, &mut q[..n + 3])?;
        self.normal_wall(
            from_slice(&q[..n + 3], n),
            p,
            outward_normal,
            out,
            &mut cancelled,
        )
    }

    fn normal_wall(
        &self,
        inside: NormalCell<'_>,
        primitive: Primitive,
        outward_normal: f64,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let n = self.species_count();
        if out.len() != n + 3 {
            return Err(Error::Invalid("planar mixture wall output axis"));
        }
        let mut normal = [0.0; MAX_COMPONENTS];
        let report = self.wall_flux(
            self.normal_reduction(inside)?,
            primitive,
            outward_normal,
            &mut normal[..n + 2],
            cancelled,
        )?;
        out.fill(0.0);
        out[n] = normal[n];
        Ok(report)
    }

    pub(super) fn admissible_normal_cell(&self, cell: NormalCell<'_>) -> bool {
        self.normal_reduction(cell)
            .is_ok_and(|q| self.admissible_cell(q))
    }

    pub(super) fn normal_reduction<'a>(&self, cell: NormalCell<'a>) -> Result<Cell<'a>, Error> {
        let (rho, _) = self.admit_densities(cell.partial_densities)?;
        if ![cell.normal_momentum, cell.transverse_momentum, cell.energy]
            .iter()
            .all(|x| x.is_finite())
        {
            return Err(Error::Invalid("non-finite planar mixture trace"));
        }
        let energy =
            cell.energy - 0.5 * cell.transverse_momentum * (cell.transverse_momentum / rho);
        if !energy.is_finite() {
            return Err(Error::Invalid(
                "unrepresentable planar mixture kinetic energy",
            ));
        }
        Ok(Cell {
            partial_densities: cell.partial_densities,
            momentum: cell.normal_momentum,
            energy,
        })
    }
}

fn from_slice(q: &[f64], n: usize) -> NormalCell<'_> {
    NormalCell {
        partial_densities: &q[..n],
        normal_momentum: q[n],
        transverse_momentum: q[n + 1],
        energy: q[n + 2],
    }
}
fn write_cell(cell: NormalCell<'_>, out: &mut [f64]) {
    let n = cell.partial_densities.len();
    out[..n].copy_from_slice(cell.partial_densities);
    out[n] = cell.normal_momentum;
    out[n + 1] = cell.transverse_momentum;
    out[n + 2] = cell.energy;
}

#[cfg(test)]
mod tests;
