//! Partition original matter without requerying density or changing cell mass.

use super::*;

/// A compact part of the original material mesh. The index maps give the
/// original source vertex/cell for every entry in `body`, in ascending order.
/// They carry nodal fields and per-cell history labels without treating those
/// fields as new material queries. Constitutive history transport remains with
/// its law; this owner only retains the original geometry and density state.
#[derive(Clone, Debug, PartialEq)]
pub struct TetrahedralDensityPart {
    /// Original material cells and sources in compact local topology.
    pub body: TetrahedralDensityBody,
    /// Original vertex index for each local vertex and its nodal fields.
    pub original_vertex_indices: Vec<usize>,
    /// Original canonical cell index for each local cell and its history.
    pub original_cell_indices: Vec<usize>,
}

/// Two nonempty material pieces containing every original cell exactly once.
/// Shared interface vertices occur in both local meshes; their copies do not
/// duplicate material mass, which remains attached to the original cells.
#[derive(Clone, Debug, PartialEq)]
pub struct TetrahedralDensitySplit {
    /// Exactly the requested original cells, in original canonical order.
    pub selected: TetrahedralDensityPart,
    /// Every original cell not selected, in original canonical order.
    pub retained: TetrahedralDensityPart,
}

fn reserved<T>(count: usize) -> Result<Vec<T>, MaterialInertiaError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| invalid("material partition allocation refused"))?;
    Ok(values)
}

fn source_vertices(
    source: &TetrahedralDensityBody,
    cells: &[usize],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<usize>, MaterialInertiaError> {
    // Each original cell contributes four references. This fixed hard bound
    // also bounds sorting between the two cancellation checkpoints.
    let mut vertices = reserved(4 * cells.len())?;
    for &cell in cells {
        poll(cancelled)?;
        vertices.extend_from_slice(&source.tetrahedra[cell]);
    }
    vertices.sort_unstable();
    vertices.dedup();
    poll(cancelled)?;
    Ok(vertices)
}

fn assemble(
    source: &TetrahedralDensityBody,
    original_cell_indices: Vec<usize>,
    original_vertex_indices: Vec<usize>,
    local_indices: &mut [usize],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<TetrahedralDensityPart, MaterialInertiaError> {
    let mut vertices = reserved(original_vertex_indices.len())?;
    for (local, &original) in original_vertex_indices.iter().enumerate() {
        poll(cancelled)?;
        local_indices[original] = local;
        vertices.push(source.vertices[original].clone());
    }
    let mut tetrahedra = reserved(original_cell_indices.len())?;
    let mut determinants = reserved(original_cell_indices.len())?;
    let mut cells = reserved(original_cell_indices.len())?;
    let mut combined = None;
    for &original in &original_cell_indices {
        poll(cancelled)?;
        // Monotone vertex renumbering preserves canonical cell order and the
        // density/determinant convention of the exact retained cell payload.
        tetrahedra.push(source.tetrahedra[original].map(|v| local_indices[v]));
        determinants.push(source.determinants[original]);
        let cell = source.cells[original].clone();
        combined = Some(combine(combined, cell.inertia)?);
        cells.push(cell);
    }
    let inertia = combined.expect("partition admission requires both cell sets nonempty");
    inertia.matrix()?;
    poll(cancelled)?;
    Ok(TetrahedralDensityPart {
        body: TetrahedralDensityBody {
            vertices,
            tetrahedra,
            determinants,
            cells,
            inertia,
        },
        original_vertex_indices,
        original_cell_indices,
    })
}

impl TetrahedralDensityBody {
    /// Split selected original cells from the retained material, with compact
    /// shared-vertex topology and explicit source maps for both parts.
    ///
    /// `selected` uses `tetrahedra()` order. Selection order is immaterial;
    /// duplicates, invalid indices and empty/full selections refuse. Each cell
    /// belongs to exactly one output. Every resolved vertex source receipt,
    /// original cell mass/density/full inertia and orientation determinant is
    /// retained unchanged. Native inertia merging supplies each new complete
    /// body. There is no density requery, cell integration or nominal-density
    /// approximation, and no source mutation on any refusal.
    ///
    /// `maximum_cells` bounds the complete original cell count. The second
    /// limit counts vertices across BOTH compact outputs, including interface
    /// copies; it is checked before cloning source payloads. Existing hard
    /// limits remain binding. Scratch and copied topology are O(vertices +
    /// cells); source-index sorting costs O(cells log cells). Cancellation is
    /// checked per cell/vertex and around each bounded sort/native merge.
    ///
    /// This partitions material geometry. It neither evolves momentum nor
    /// admits separation of a bond, fracture energy, phase change, or elastic
    /// and plastic history deletion. A physical transition must supply those
    /// laws and transfer the original motion and fields through the maps.
    pub fn partition_cells(
        &self,
        selected: &[usize],
        maximum_cells: usize,
        maximum_output_vertices: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TetrahedralDensitySplit, MaterialInertiaError> {
        poll(&mut cancelled)?;
        let n = self.tetrahedra.len();
        if maximum_cells == 0
            || maximum_cells > MAX_BODY_REGIONS
            || n > maximum_cells
            || maximum_output_vertices == 0
            || maximum_output_vertices > MAX_DENSITY_VERTICES
            || selected.is_empty()
            || selected.len() >= n
        {
            return Err(invalid(
                "bounded source and two nonempty material cell sets required",
            ));
        }
        let mut chosen = reserved(n)?;
        chosen.resize(n, false);
        for &index in selected {
            poll(&mut cancelled)?;
            if index >= n || chosen[index] {
                return Err(invalid(
                    "selected material cells must be valid and distinct",
                ));
            }
            chosen[index] = true;
        }
        let mut selected_cells = reserved(selected.len())?;
        let mut retained_cells = reserved(n - selected.len())?;
        for (cell, &is_selected) in chosen.iter().enumerate() {
            poll(&mut cancelled)?;
            if is_selected {
                selected_cells.push(cell);
            } else {
                retained_cells.push(cell);
            }
        }
        let selected_vertices = source_vertices(self, &selected_cells, &mut cancelled)?;
        let retained_vertices = source_vertices(self, &retained_cells, &mut cancelled)?;
        if selected_vertices.len() + retained_vertices.len() > maximum_output_vertices {
            return Err(invalid(
                "material partition exceeds its combined output-vertex ceiling",
            ));
        }
        let mut local_indices = reserved(self.vertices.len())?;
        local_indices.resize(self.vertices.len(), usize::MAX);
        let selected = assemble(
            self,
            selected_cells,
            selected_vertices,
            &mut local_indices,
            &mut cancelled,
        )?;
        let retained = assemble(
            self,
            retained_cells,
            retained_vertices,
            &mut local_indices,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(TetrahedralDensitySplit { selected, retained })
    }
}
