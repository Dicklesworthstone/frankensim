// A visualization boundary, not a mesh repairer, physical solver or certificate.
export const MESH_SCHEMA = 'frankensim.observatory/mesh-v1';
export const MESH_LIMITS = Object.freeze({bytes: 32 * 1024 * 1024, points: 50000,
  cells: 50000, triangles: 100000, frames: 512, fields: 24, values: 2000000});
export const fail = message => {throw new Error(message);};
export function scalar(value, label = 'Value') {
  if (typeof value !== 'number' || !Number.isFinite(value) || Math.abs(value) > 1e100) fail(`${label} must be finite and at most 1e100 in magnitude.`);
  return value;
}
export function label(value, name = 'Label', limit = 200) {
  if (typeof value !== 'string' || !value.trim() || value.length > limit) fail(`Invalid ${name}.`);
  return value;
}
function array(value, name, min, max) {
  if (!Array.isArray(value) || value.length < min || value.length > max) fail(`${name}: expected ${min}–${max} entries.`);
  return value;
}
function integer(value, name, min, max) {
  if (!Number.isSafeInteger(value) || value < min || value > max) fail(`Invalid ${name}.`);
  return value;
}
export function boundsOf(points) {
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  for (const p of points) for (let k = 0; k < 3; k++) {min[k] = Math.min(min[k], p[k]); max[k] = Math.max(max[k], p[k]);}
  const center = min.map((x, k) => x + (max[k] - x) / 2);
  return {min, max, center, span: Math.max(...max.map((x, k) => x - min[k]))};
}
export function validateMeshSnapshot(raw) {
  if (!raw || raw.schema !== MESH_SCHEMA) fail(`Expected ${MESH_SCHEMA}.`);
  const source = raw.source;
  if (!source) fail('Source description is required.');
  const result = {schema: MESH_SCHEMA, title: label(raw.title, 'title'), source: {
    model: label(source.model, 'model'), material: label(source.material, 'material'),
    notes: array(source.notes, 'Source notes', 1, 16).map(n => label(n, 'source note', 1500))}};
  const m = raw.mesh; if (!m) fail('Mesh is required.');
  const positions_m = array(m.positions_m, 'Positions', 3, MESH_LIMITS.points).map(p =>
    array(p, 'XYZ coordinate', 3, 3).map(v => {
      scalar(v, 'Coordinate'); if (Math.abs(v) > 1e30) fail('Coordinate exceeds display range (1e30 m).'); return v;
    }));
  const bounds = boundsOf(positions_m);
  if (!(bounds.span > 0)) fail('Mesh has zero extent.');
  const cell_count = integer(m.cell_count, 'source cell count', 1, MESH_LIMITS.cells);
  const used = new Set();
  const triangles = array(m.triangles, 'Surface triangles', 1, MESH_LIMITS.triangles).map(t => {
    const tri = array(t, 'Triangle', 3, 3).map(i => integer(i, 'point index', 0, positions_m.length - 1));
    const key = [...tri].sort((a,b) => a-b).join(',');
    if (new Set(tri).size !== 3 || used.has(key)) fail('Repeated or duplicate surface triangle.');
    used.add(key);
    const a = positions_m[tri[0]], b = positions_m[tri[1]], c = positions_m[tri[2]];
    const u = b.map((x,k) => (x-a[k])/bounds.span), v = c.map((x,k) => (x-a[k])/bounds.span);
    if (Math.hypot(u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]) === 0) fail('Degenerate surface triangle.');
    return tri;
  });
  const triangle_cells = array(m.triangle_cells, 'Surface/source mapping', triangles.length, triangles.length)
    .map(i => integer(i, 'source cell index', 0, cell_count-1));
  const time_s = array(m.time_s, 'Recorded frames', 1, MESH_LIMITS.frames).map(v => scalar(v, 'Frame time'));
  if (time_s[0] < 0 || time_s.some((v,i) => i > 0 && v <= time_s[i-1])) fail('Recorded times must be nonnegative and strictly increasing.');
  let count = positions_m.length * 3 + triangles.length * 4;
  let source_cells, source_types;
  if (m.source_cells != null || m.source_types != null) {
    source_types = array(m.source_types, 'Original cell types', cell_count, cell_count).map(t => {
      if (![5,9,10,12].includes(t)) fail('Unsupported original cell type.'); return t;
    });
    source_cells = array(m.source_cells, 'Original cell connectivity', cell_count, cell_count).map((cell,i) => {
      const size = {5:3,9:4,10:4,12:8}[source_types[i]];
      const ids = array(cell, 'Original cell', size, size).map(v => integer(v, 'original cell point', 0, positions_m.length-1));
      if (new Set(ids).size !== size) fail('Repeated point in original cell.'); count += size+1; return ids;
    });
    triangles.forEach((tri,i) => {
      if (tri.some(p => !source_cells[triangle_cells[i]].includes(p))) fail('Surface face does not belong to its declared original cell.');
    });
  }
  if (count > MESH_LIMITS.values) fail('Mesh exceeds the 2,000,000-value working-set budget.');
  const names = new Set();
  const fields = array(m.fields, 'Fields', 0, MESH_LIMITS.fields).map(f => {
    const name = label(f.label, 'field label', 100), unit = label(f.unit, 'field unit', 32);
    if (!['point', 'cell'].includes(f.association) || ![1,3].includes(f.components)) fail('Fields require point/cell association and 1 or 3 components.');
    const key = `${f.association}:${name}`; if (names.has(key)) fail('Duplicate field name within one association.'); names.add(key);
    const n = (f.association === 'point' ? positions_m.length : cell_count) * f.components;
    count += n * time_s.length; if (count > MESH_LIMITS.values) fail('Mesh exceeds the 2,000,000-value working-set budget.');
    const values = array(f.values, 'Field frames', time_s.length, time_s.length)
      .map(row => array(row, 'Field values', n, n).map(v => scalar(v, name)));
    return {label: name, unit, association: f.association, components: f.components, values};
  });
  const mesh = {positions_m, triangles, triangle_cells, cell_count, time_s, fields};
  if (source_cells) {mesh.source_cells = source_cells; mesh.source_types = source_types;}
  if (m.displacement_field != null) {
    const index = integer(m.displacement_field, 'displacement field index', 0, fields.length-1), f = fields[index];
    if (f.association !== 'point' || f.components !== 3 || f.unit !== 'm') fail('Displacement must explicitly select a three-component point field in metres.');
    mesh.displacement_field = index;
  }
  result.mesh = mesh; return result;
}
// A recorded-frame hold, not interpolation or a new integration scheme.
export function frameAt(times, time) {
  let lo = 0, hi = times.length;
  while (lo < hi) {const mid = (lo+hi) >>> 1; if (times[mid] <= time) lo = mid+1; else hi = mid;}
  return Math.max(0, lo-1);
}
export function fieldValue(field, frame, entity, component = 'magnitude') {
  const row = field.values[frame], i = entity * field.components;
  return field.components === 1 ? row[i] : component === 'magnitude'
    ? Math.hypot(row[i], row[i+1], row[i+2]) : row[i + Number(component)];
}
export function fieldRange(field, component = 'magnitude') {
  let min = Infinity, max = -Infinity;
  for (let frame = 0; frame < field.values.length; frame++) {
    for (let i = 0; i < field.values[frame].length/field.components; i++) {
      const v = fieldValue(field, frame, i, component); min = Math.min(min,v); max = Math.max(max,v);
    }
  }
  return [min,max];
}
export function probe(mesh, fieldIndex, entity, component = 'magnitude') {
  const f = mesh.fields[fieldIndex]; if (!f) return null;
  integer(entity, 'probe entity', 0, (f.association === 'point' ? mesh.positions_m.length : mesh.cell_count)-1);
  return {time_s: mesh.time_s, series: [{id:'probe', label:`${f.label} · ${f.association} ${entity}${f.components===3 ? ' · '+(component==='magnitude'?'magnitude':'xyz'[Number(component)]) : ''}`,
    unit: f.unit, values: mesh.time_s.map((_,i) => fieldValue(f,i,entity,component))}]};
}
// No topology resampling: a sequence must retain exactly the same entity identities.
export function joinFrames(snapshots, times) {
  if (snapshots.length !== times.length || snapshots.length < 1 || snapshots.length > MESH_LIMITS.frames) fail('One explicit time is required per snapshot.');
  const base = validateMeshSnapshot(snapshots[0]);
  const topology = m => JSON.stringify([m.positions_m,m.triangles,m.triangle_cells,m.cell_count,m.displacement_field,m.source_cells,m.source_types]);
  const shape = m => JSON.stringify(m.fields.map(f => [f.label,f.unit,f.association,f.components]));
  const signature = topology(base.mesh), fields = shape(base.mesh);
  if (base.mesh.fields.some(f => f.association === 'cell') && !base.mesh.source_cells) fail('Cell-valued sequences require original source-cell connectivity.');
  for (let i=0;i<snapshots.length;i++) {
    const s = i===0 ? base : validateMeshSnapshot(snapshots[i]);
    if (s.mesh.time_s.length!==1 || topology(s.mesh)!==signature || shape(s.mesh)!==fields
      || JSON.stringify(s.source)!==JSON.stringify(base.source)) fail('Sequence changed geometry, field layout or source; import these as separate runs.');
  }
  base.mesh.time_s = [...times];
  base.mesh.fields = base.mesh.fields.map((f,j) => ({...f,values:snapshots.map(s => s.mesh.fields[j].values[0])}));
  return validateMeshSnapshot(base);
}
