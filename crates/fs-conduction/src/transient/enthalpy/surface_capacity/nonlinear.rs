//! Original total film inventories distributed over their one boundary face.

use super::*;
use fs_material::thermal::transport::{
    ThermalTransportError,
    storage::{PreparedHeatStorage, TabulatedHeatStorage},
};

#[derive(Debug, Clone, Copy)]
enum CaloricStorage<'s> {
    Prepared(&'s PreparedHeatStorage),
    Tabulated(&'s TabulatedHeatStorage),
}

/// Nonlinear sensible storage of a uniform, perfectly bonded boundary film.
///
/// The original law describes the TOTAL inventory of this one triangle, in
/// J/K, not an areal capacity. Its fixed amounts are distributed by the three
/// equal P1 row moments, each one third of that inventory. At nonuniform host
/// temperatures this is nodal lumping, not exact nonlinear surface quadrature.
///
/// Each face retains its original borrowed caloric object. Reusing that object
/// on a second face explicitly declares a second equal inventory. No material,
/// mass, reference area, heat-capacity constraint or source authority is inferred.
/// A separately meshed coating must not be added again. The film shares the
/// host temperature and has no independent phase state or latent heat.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceHeatStorage<'s> {
    boundary_face: usize,
    law: CaloricStorage<'s>,
}

impl<'s> SurfaceHeatStorage<'s> {
    /// Bind the original source-qualified total inventory of this face.
    #[must_use]
    pub const fn prepared(boundary_face: usize, storage: &'s PreparedHeatStorage) -> Self {
        Self {
            boundary_face,
            law: CaloricStorage::Prepared(storage),
        }
    }

    /// Bind a caller-declared total C(T) table, without material-source authority.
    #[must_use]
    pub const fn tabulated(boundary_face: usize, storage: &'s TabulatedHeatStorage) -> Self {
        Self {
            boundary_face,
            law: CaloricStorage::Tabulated(storage),
        }
    }

    /// Slot in the original mesh's boundary-face table.
    #[must_use]
    pub const fn boundary_face(self) -> usize {
        self.boundary_face
    }

    /// The exact source-qualified inventory, when explicitly supplied.
    #[must_use]
    pub const fn prepared_storage(self) -> Option<&'s PreparedHeatStorage> {
        match self.law {
            CaloricStorage::Prepared(storage) => Some(storage),
            CaloricStorage::Tabulated(_) => None,
        }
    }

    /// The exact declared table, when explicitly supplied.
    #[must_use]
    pub const fn tabulated_storage(self) -> Option<&'s TabulatedHeatStorage> {
        match self.law {
            CaloricStorage::Tabulated(storage) => Some(storage),
            CaloricStorage::Prepared(_) => None,
        }
    }

    fn temperature_span(self) -> TemperatureSpan {
        let [low, high] = match self.law {
            CaloricStorage::Prepared(storage) => storage.temperature_span_k(),
            CaloricStorage::Tabulated(storage) => storage.temperature_span_k(),
        };
        TemperatureSpan::Sampled { low, high }
    }

    pub(super) fn heat_capacity_j_k(self, temperature: f64) -> Result<f64, ThermalTransportError> {
        match self.law {
            CaloricStorage::Prepared(storage) => storage.heat_capacity_j_k(temperature),
            CaloricStorage::Tabulated(storage) => storage.heat_capacity_j_k(temperature),
        }
    }

    pub(super) fn energy_change_j(self, old: f64, new: f64) -> Result<f64, ThermalTransportError> {
        match self.law {
            CaloricStorage::Prepared(storage) => storage.energy_change_j(old, new),
            CaloricStorage::Tabulated(storage) => storage.energy_change_j(old, new),
        }
    }
}

#[derive(Debug, Default)]
pub(super) struct NonlinearSurfaceStorage<'s> {
    faces: Vec<SurfaceHeatStorage<'s>>,
    // Compiled sparse incidence. Constitutive queries never scan unrelated
    // faces, allocate, clone tables, or resolve source claims again.
    offsets: Vec<usize>,
    indices: Vec<usize>,
}

impl<'s> NonlinearSurfaceStorage<'s> {
    pub(super) fn at_vertex(&self, vertex: usize) -> impl Iterator<Item = &SurfaceHeatStorage<'s>> {
        self.indices[self.offsets[vertex]..self.offsets[vertex + 1]]
            .iter()
            .map(|&index| &self.faces[index])
    }

    fn rebuild(&mut self, cx: &Cx<'_>, mesh: &ConductionMesh) -> Result<(), EnthalpyError> {
        let count = self
            .faces
            .len()
            .checked_mul(3)
            .filter(|&count| count <= isize::MAX as usize / std::mem::size_of::<usize>())
            .ok_or(EnthalpyError::InvalidInput(
                "surface-storage incidence size overflow",
            ))?;
        self.offsets = vec![0; mesh.vertex_count() + 1];
        for (index, face) in self.faces.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            for vertex in mesh.boundary()[face.boundary_face].vertices {
                self.offsets[vertex as usize + 1] += 1;
            }
        }
        for vertex in 0..mesh.vertex_count() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            self.offsets[vertex + 1] += self.offsets[vertex];
        }
        let mut cursor = self.offsets[..mesh.vertex_count()].to_vec();
        self.indices = vec![0; count];
        for (index, face) in self.faces.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            for vertex in mesh.boundary()[face.boundary_face].vertices {
                let slot = &mut cursor[vertex as usize];
                self.indices[*slot] = index;
                *slot += 1;
            }
        }
        poll(cx, count)
    }
}

impl<'m, 'c> EnthalpyBackwardEuler<'m, 'c> {
    /// Add nonlinear sensible films to the original host's primary enthalpy.
    ///
    /// Each supplied law owns one face's TOTAL inventory. A receiving row adds
    /// one third of its original integral C(T)dT between old/new temperatures;
    /// its tangent adds C(T)/3 times the host's dT/dh. Thus old, stage and endpoint
    /// slopes are evaluated at their own temperatures, and a host latent plateau
    /// stores no sensible film heat. Tables are borrowed, never sampled/frozen.
    ///
    /// BE, SDIRK2, adaptivity, contact/radiation, bath reactions and adjoints use
    /// the same storage helpers. Only the immutable data path is prepared here;
    /// no source query or allocation occurs in a film's caloric evaluation.
    /// Coating-parameter derivatives and independent film phase changes are
    /// outside this lane. The caller admits the original storage work constraint.
    ///
    /// Faces may appear once per call; further calls explicitly add further
    /// layers. Invalid faces, duplicate faces, disjoint original source support
    /// or cancellation publish no partly prepared model. Both history and trial
    /// temperatures must remain inside every original participating film's support.
    pub fn with_surface_heat_storage(
        mut self,
        cx: &Cx<'_>,
        surfaces: &[SurfaceHeatStorage<'c>],
    ) -> Result<Self, EnthalpyError> {
        poll(cx, 0)?;
        let vertex_count = self.masses.len();
        let mut seen = BTreeSet::new();
        for (index, &surface) in surfaces.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            let face = self.mesh.boundary().get(surface.boundary_face).ok_or(
                EnthalpyError::InvalidInput(
                    "surface storage face is outside the original boundary table",
                ),
            )?;
            if !seen.insert(surface.boundary_face) {
                return Err(EnthalpyError::InvalidInput(
                    "surface storage face is declared more than once",
                ));
            }
            let vertices = face.vertices.map(|vertex| vertex as usize);
            let spans = vertices.map(|vertex| {
                let knots = self.curve_for_vertex(vertex).knots();
                surface
                    .temperature_span()
                    .intersect(TemperatureSpan::Sampled {
                        low: knots[0].temperature_k,
                        high: knots[knots.len() - 1].temperature_k,
                    })
            });
            let prepared = self
                .surface_capacity
                .get_or_insert_with(|| SurfaceStorage::new(vertex_count));
            for (&vertex, span) in vertices.iter().zip(spans) {
                let span = prepared.temperature_spans[vertex].intersect(span);
                if matches!(span, TemperatureSpan::Sampled { low, high } if low > high) {
                    return Err(EnthalpyError::InvalidInput(
                        "host chart and film storage sources have disjoint temperature support",
                    ));
                }
                prepared.temperature_spans[vertex] = span;
            }
            prepared
                .nonlinear
                .get_or_insert_with(NonlinearSurfaceStorage::default)
                .faces
                .push(surface);
        }
        if !surfaces.is_empty() {
            self.surface_capacity
                .as_mut()
                .expect("active films prepared")
                .nonlinear
                .as_mut()
                .expect("nonlinear films prepared")
                .rebuild(cx, self.mesh)?;
        }
        poll(cx, surfaces.len())?;
        Ok(self)
    }

    /// Original nonlinear face inventories, preserving source/table identity.
    #[must_use]
    pub fn surface_heat_storage(&self) -> &[SurfaceHeatStorage<'c>] {
        self.surface_capacity
            .as_ref()
            .and_then(|film| film.nonlinear.as_ref())
            .map_or(&[], |nonlinear| nonlinear.faces.as_slice())
    }
}
