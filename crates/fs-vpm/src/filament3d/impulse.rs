//! Geometric angular impulse of the same closed filament rings as source assembly.

use super::{FilamentWake, Refusal, ring_vector_area};
use fs_exec::reduce::Compensated;

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Unit-circulation `-1/2 ∮ |x-origin|² dx` on one oriented quad.
/// A triangle contributes its centroid arm crossed with its vector area.
/// Four triangles around a nearby anchor avoid subtracting squared world
/// coordinates; opposite-edge pairing preserves reflected-sheet symmetry.
/// The anchor transport uses the existing diagonal-based vector-area owner.
pub(crate) fn ring_angular_impulse(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    origin: [f64; 3],
) -> [f64; 3] {
    let anchor: [f64; 3] =
        core::array::from_fn(|i| (0.25 * a[i] + 0.25 * c[i]) + (0.25 * b[i] + 0.25 * d[i]));
    let offsets = [a, b, c, d].map(|p| core::array::from_fn(|i| p[i] - anchor[i]));
    let terms: [[f64; 3]; 4] = core::array::from_fn(|edge| {
        let p = offsets[edge];
        let q = offsets[(edge + 1) % 4];
        let centroid = core::array::from_fn(|i| (p[i] + q[i]) / 3.0);
        cross(centroid, cross(p, q).map(|v| 0.5 * v))
    });
    let transport = cross(
        core::array::from_fn(|i| anchor[i] - origin[i]),
        ring_vector_area(a, b, c, d),
    );
    core::array::from_fn(|i| {
        ((terms[0][i] + terms[2][i]) + (terms[1][i] + terms[3][i])) + transport[i]
    })
}

impl FilamentWake {
    /// Geometric hydrodynamic angular impulse about `origin`, divided by
    /// uniform fluid density [m⁵/s]. Each closed filament ring contributes
    /// `-Γ/2 ∮ |x-origin|² dx`, with the same orientation as source assembly.
    /// Moving the origin by `a` gives `J(origin+a) = J(origin) - a × I`.
    /// Ring moments are evaluated locally before this origin transport, so a
    /// remote world origin does not require cancelling squared coordinates.
    ///
    /// This is the localized-vorticity angular impulse (Carpenter & Guha,
    /// Phys. Rev. Fluids 5, 064702, 2020, doi:10.1103/PhysRevFluids.5.064702),
    /// not the angular momentum inside an arbitrary finite fluid control volume.
    /// Core smoothing does not change this geometric filament diagnostic.
    ///
    /// # Errors
    /// `impulse-interface-attached` for an open retained near/far partition;
    /// only the combined closed sources have this origin-transport law.
    /// Invalid geometry, a nonfinite origin or an unrepresentable moment refuses.
    pub fn hydrodynamic_angular_impulse_per_density(
        &self,
        origin: [f64; 3],
    ) -> Result<[f64; 3], Refusal> {
        if self.downstream_gamma.is_some() {
            return Err(Refusal {
                code: "impulse-interface-attached",
                message: "retained near sources alone do not form a complete closed wake".into(),
                ranked_repairs: vec!["measure the complete wake before near/far partition".into()],
            });
        }
        self.admit_geometry()?;
        if origin.iter().any(|v| !v.is_finite()) {
            return Err(Refusal {
                code: "impulse-origin-invalid",
                message: "angular impulse needs a finite reference origin".into(),
                ranked_repairs: vec!["choose a finite reference origin".into()],
            });
        }
        let mut result = [Compensated::zero(); 3];
        for (r, row) in self.rows.iter().enumerate() {
            let next = self
                .rows
                .get(r + 1)
                .map_or(&self.line_nodes, |row| &row.nodes);
            for s in 0..self.cert.n_stations {
                let moment = ring_angular_impulse(
                    next[s],
                    next[s + 1],
                    row.nodes[s + 1],
                    row.nodes[s],
                    origin,
                );
                for i in 0..3 {
                    result[i] = result[i].accumulate(row.gamma[s] * moment[i]);
                }
            }
        }
        let result = result.map(Compensated::value);
        if result.iter().any(|v| !v.is_finite()) {
            return Err(Refusal {
                code: "impulse-unrepresentable",
                message: "circulation times angular moment exceeds the numerical range".into(),
                ranked_repairs: vec![
                    "rescale geometry and circulation or use a nearby origin".into(),
                ],
            });
        }
        Ok(result)
    }
}
