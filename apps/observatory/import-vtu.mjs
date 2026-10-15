// Narrow offline reader for fs-viz::VtuWriter::write_ascii and compatible VTU.
// No external entities, remote references, binary payloads or XML evaluation.
import {MESH_SCHEMA, MESH_LIMITS, fail, scalar, validateMeshSnapshot} from './mesh-data.mjs';
const numeric = /^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$/;
function decode(s) {
  return s.replace(/&([^;]*);/g, (_, code) => {
    const named = {amp:'&',lt:'<',gt:'>',quot:'"',apos:"'"};
    if (Object.hasOwn(named,code)) return named[code];
    if (/^#(?:[0-9]+|x[0-9a-f]+)$/i.test(code)) {
      const n = code[1].toLowerCase()==='x' ? parseInt(code.slice(2),16) : Number(code.slice(1));
      if (n>0 && n<=0x10ffff && !(n>=0xd800 && n<=0xdfff)) return String.fromCodePoint(n);
    }
    fail('Unsupported XML entity.');
  });
}
function xml(source) {
  const root = {tag:'document',children:[],text:'',attrs:{}}, stack = [root];
  const lexer = /<!--[\s\S]*?-->|<\?xml\s[^?]*\?>|<(?:"[^"]*"|'[^']*'|[^'">])*>|[^<]+/g;
  let end=0, count=0;
  for (const match of source.matchAll(lexer)) {
    if (match.index!==end || ++count>1000) fail('Malformed or overly complex VTU XML.');
    end=match.index+match[0].length; const token=match[0], parent=stack.at(-1);
    if (token.startsWith('<!--')) continue;
    if (token.startsWith('<?xml')) {if (parent!==root || root.children.length) fail('Unexpected XML declaration.'); continue;}
    if (!token.startsWith('<')) {parent.text+=token; continue;}
    const closing=token.match(/^<\/([A-Za-z][\w]*)\s*>$/);
    if (closing) {if (stack.length<2 || parent.tag!==closing[1]) fail('Mismatched XML element.'); stack.pop(); continue;}
    const opening=token.match(/^<([A-Za-z][\w]*)([\s\S]*?)(\/?)>$/);
    if (!opening || stack.length>7) fail('Unsupported VTU XML declaration or nesting.');
    const node={tag:opening[1],children:[],text:'',attrs:{}}; let attributes=opening[2];
    while (attributes.trim()) {
      const attr=attributes.match(/^\s+([A-Za-z_][\w.-]*)\s*=\s*(?:"([^"<]*)"|'([^'<]*)')/);
      if (!attr || Object.hasOwn(node.attrs,attr[1])) fail('Malformed or duplicate XML attribute.');
      const value=attr[2]??attr[3];
      if (/&(?!(?:amp|lt|gt|quot|apos|#\d+|#x[0-9a-fA-F]+);)/.test(value)) fail('Malformed XML entity.');
      node.attrs[attr[1]]=decode(value); attributes=attributes.slice(attr[0].length);
    }
    parent.children.push(node); if (!opening[3]) stack.push(node);
  }
  if (end!==source.length || stack.length!==1 || root.children.length!==1 || root.text.trim()) fail('Incomplete VTU XML.');
  return root.children[0];
}
function only(node, allowed) {
  if (node.text.trim() || node.children.some(c=>!allowed.includes(c.tag))) fail(`Unsupported content inside ${node.tag}.`);
}
function one(node, tag, optional=false) {
  const all=node.children.filter(c=>c.tag===tag);
  if (all.length!==1 && !(optional && all.length===0)) fail(`Expected one ${tag}.`);
  return all[0];
}
function uint(text, name, max) {
  if (!/^\d+$/.test(text||'') || !Number.isSafeInteger(Number(text)) || Number(text)>max) fail(`Invalid ${name}.`);
  return Number(text);
}
function data(node, count, {integer=false, components=null}={}) {
  if (!node || node.tag!=='DataArray' || node.children.length || node.attrs.format!=='ascii') fail('Only inline ASCII DataArray values are supported.');
  if (!['Float64','Float32','Int32','Int64','UInt8','UInt32','UInt64'].includes(node.attrs.type)) fail('Unsupported VTU numeric type.');
  if (components!==null && uint(node.attrs.NumberOfComponents||'1','component count',9)!==components) fail('Incorrect VTU component count.');
  const text=node.text.trim(), tokens=text ? text.split(/\s+/) : [];
  if (tokens.length!==count) fail(`DataArray ${node.attrs.Name||'Points'}: expected ${count} values, found ${tokens.length}.`);
  return tokens.map(t=>{
    if (!numeric.test(t)) fail('Invalid numeric token in VTU.');
    const x=scalar(Number(t),'VTU value');
    if ((integer || !node.attrs.type.startsWith('Float')) && !Number.isSafeInteger(x)) fail('VTU integer exceeds JavaScript exact-integer range.');
    return x;
  });
}
const sub=(a,b)=>a.map((v,i)=>v-b[i]);
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const dot=(a,b)=>a.reduce((sum,v,i)=>sum+v*b[i],0);
export function exterior(points,cells,types) {
  const faces=new Map(), sourceCells=new Set(); let volume=null;
  const centroid=ids=>[0,1,2].map(k=>ids.reduce((sum,id)=>sum+points[id][k]/ids.length,0));
  for (let cell=0;cell<cells.length;cell++) {
    const ids=cells[cell], type=types[cell];
    const size={5:3,9:4,10:4,12:8}[type];
    if (!size || ids.length!==size) fail(`Cell ${cell}: only triangles, quads, tetrahedra and hexahedra are supported.`);
    if (ids.some(i=>!Number.isSafeInteger(i)||i<0||i>=points.length) || new Set(ids).size!==ids.length) fail(`Cell ${cell}: invalid or repeated point index.`);
    const solid=type===10||type===12;
    if (volume!==null && volume!==solid) fail('Mixed volume/surface cells are ambiguous; export either the volume or its surface.');
    volume=solid;
    const identity=[...ids].sort((a,b)=>a-b).join(','); if (sourceCells.has(identity)) fail('Duplicate source cell.'); sourceCells.add(identity);
    const local=type===10 ? [[0,2,1],[0,1,3],[1,2,3],[2,0,3]] : type===12
      ? [[0,3,2,1],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]] : [ids.map((_,i)=>i)];
    for (const face of local) {
      const f=face.map(i=>ids[i]), a=points[f[0]], u=sub(points[f[1]],a), v=sub(points[f[2]],a);
      const scale=Math.max(...u.map(Math.abs),...v.map(Math.abs));
      if (!(scale>0)) fail('Degenerate cell face.');
      const normal=cross(u.map(x=>x/scale),v.map(x=>x/scale));
      if (Math.hypot(...normal)===0) fail('Degenerate cell face.');
      if (f.length===4) {
        const w=sub(points[f[3]],a).map(x=>x/scale);
        if (Math.abs(dot(normal,w))>1e-8*Math.hypot(...normal)*Math.max(1,Math.hypot(...w))) fail('Nonplanar quad face: triangulate explicitly before importing.');
        for(let i=0;i<4;i++) {
          const p=points[f[i]], b=sub(points[f[(i+1)%4]],p).map(x=>x/scale), c=sub(points[f[(i+2)%4]],points[f[(i+1)%4]]).map(x=>x/scale);
          if (!(dot(cross(b,c),normal)>0)) fail('Nonconvex or degenerate quad face.');
        }
      }
      if (solid) {
        const side=dot(normal,sub(centroid(ids),a).map(x=>x/scale));
        if (side===0) fail('Degenerate volume cell.');
        if (side>0) f.reverse();
        // This is a convex-cell display admission, not global mesh certification.
        const n=side>0?normal.map(x=>-x):normal;
        for(const id of ids) if (!f.includes(id) && !(dot(n,sub(points[id],a).map(x=>x/scale))<0)) fail('Nonconvex volume cell.');
      }
      const key=[...f].sort((a,b)=>a-b).join(','), old=faces.get(key);
      if (old) {
        if (!solid || ++old.count>2) fail('Duplicate surface or nonmanifold volume face.');
        const index=old.ids.indexOf(f[0]);
        if (old.ids[(index+1)%f.length]!==f.at(-1)) fail('Overlapping or inconsistently joined volume cells.');
      } else faces.set(key,{ids:f,cell,count:1});
    }
  }
  const triangles=[], triangle_cells=[];
  for (const face of faces.values()) if (face.count===1) {
    for (let j=1;j<face.ids.length-1;j++) {
      if (triangles.length>=MESH_LIMITS.triangles) fail('Exterior exceeds the triangle display budget.');
      triangles.push([face.ids[0],face.ids[j],face.ids[j+1]]); triangle_cells.push(face.cell);
    }
  }
  return {triangles,triangle_cells};
}
export function parseVTU(text, {lengthUnit, time=null, title='Imported VTU fields'}={}) {
  const scales={m:1,mm:.001,cm:.01};
  if (!Object.hasOwn(scales,lengthUnit)) fail('Declare the coordinate unit explicitly: m, cm or mm.');
  if (typeof text!=='string' || new TextEncoder().encode(text).length>MESH_LIMITS.bytes) fail('VTU exceeds 32 MiB.');
  if (/<!DOCTYPE|<!ENTITY|<!\[CDATA\[/i.test(text)) fail('DTD, entities and CDATA are not supported.');
  const vtk=xml(text.replace(/^\uFEFF/,''));
  if (vtk.tag!=='VTKFile'||vtk.attrs.type!=='UnstructuredGrid'||vtk.attrs.compressor) fail('Expected uncompressed UnstructuredGrid VTU.');
  only(vtk,['UnstructuredGrid']); const grid=one(vtk,'UnstructuredGrid'); only(grid,['Piece']); const piece=one(grid,'Piece');
  only(piece,['Points','Cells','PointData','CellData']);
  const np=uint(piece.attrs.NumberOfPoints,'point count',MESH_LIMITS.points), nc=uint(piece.attrs.NumberOfCells,'cell count',MESH_LIMITS.cells);
  if (np<3||nc<1) fail('VTU has no displayable mesh.');
  const p=one(piece,'Points'); only(p,['DataArray']);
  const xyz=data(one(p,'DataArray'),np*3,{components:3});
  const positions_m=Array.from({length:np},(_,i)=>xyz.slice(i*3,i*3+3).map(v=>v*scales[lengthUnit]));
  const c=one(piece,'Cells'); only(c,['DataArray']); if(c.children.length!==3) fail('Expected connectivity, offsets and types.');
  const named=name=>{const list=c.children.filter(a=>a.attrs.Name===name); if(list.length!==1) fail(`Missing or duplicate ${name}.`); return list[0];};
  const offsets=data(named('offsets'),nc,{integer:true,components:1});
  if(offsets.some((v,i)=>v<1 || v>(i+1)*8 || (i>0 && v<=offsets[i-1]))) fail('Invalid cell offsets.');
  const connectivity=data(named('connectivity'),offsets.at(-1),{integer:true,components:1}), types=data(named('types'),nc,{integer:true,components:1});
  const cells=offsets.map((end,i)=>connectivity.slice(i?offsets[i-1]:0,end));
  const boundary=exterior(positions_m,cells,types), fields=[];
  let values=np*3+boundary.triangles.length*4;
  for (const [tag,association,entities] of [['PointData','point',np],['CellData','cell',nc]]) {
    const section=one(piece,tag,true); if(!section) continue; only(section,['DataArray']);
    for(const arr of section.children) {
      if(fields.length>=MESH_LIMITS.fields) fail('VTU exceeds the field-count budget.');
      const components=uint(arr.attrs.NumberOfComponents||'1','components',3);
      if(![1,3].includes(components)) fail('Only scalar and three-component vector fields are supported.');
      values+=entities*components; if(values>MESH_LIMITS.values) fail('VTU exceeds the working-set budget.');
      fields.push({label:arr.attrs.Name,unit:arr.attrs.Unit||'unspecified',association,components,values:[data(arr,entities*components)]});
    }
  }
  return validateMeshSnapshot({schema:MESH_SCHEMA,title,source:{model:'Imported ASCII VTU / fs-viz compatible',
    material:'Field values supplied by the imported file; source is not authenticated',notes:[
      `Coordinates explicitly interpreted as ${lengthUnit} and converted to metres. Field Unit attributes are retained without conversion.`,
      'Volume cells are shown by their exterior faces; interior cells and their original field values remain available to the probe. No cutting plane or volume rendering.',
      'Quads use a fixed diagonal for display. Point fields are linearly interpolated only for coloring; cell fields are never averaged across material interfaces.',
      time===null?'Static snapshot: zero is a viewer reference, not an asserted solver time.':`Recorded snapshot time explicitly supplied by the caller: ${time} s.`]},
    mesh:{positions_m,...boundary,cell_count:nc,time_s:[time===null?0:scalar(time,'Snapshot time')],fields}});
}
