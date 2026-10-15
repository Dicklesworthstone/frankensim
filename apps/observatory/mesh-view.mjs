// WebGL display of recorded surface/field data. CPU doubles remain authoritative
// for picking/probes; normalized float buffers are rendering intermediates only.
import {boundsOf,fieldValue,fieldRange} from './mesh-data.mjs';
import {rasterize} from './mesh-raster.mjs';
const clamp=(v,a,b)=>Math.max(a,Math.min(b,v));
export function projectPoint(p,{yaw,pitch,zoom,aspect,pan=[0,0]}) {
  const x=p[0]*Math.cos(yaw)+p[2]*Math.sin(yaw), z=-p[0]*Math.sin(yaw)+p[2]*Math.cos(yaw);
  return [x*zoom/aspect+pan[0],(p[1]*Math.cos(pitch)-z*Math.sin(pitch))*zoom+pan[1],
    p[1]*Math.sin(pitch)+z*Math.cos(pitch)];
}
export function pickTriangle(points,triangles,x,y) {
  let hit=null, nearest=-Infinity;
  triangles.forEach((ids,triangle)=>{
    const [a,b,c]=ids.map(i=>points[i]);
    const det=(b[1]-c[1])*(a[0]-c[0])+(c[0]-b[0])*(a[1]-c[1]); if(Math.abs(det)<1e-15)return;
    const u=((b[1]-c[1])*(x-c[0])+(c[0]-b[0])*(y-c[1]))/det;
    const v=((c[1]-a[1])*(x-c[0])+(a[0]-c[0])*(y-c[1]))/det, w=1-u-v;
    if(Math.min(u,v,w)<-1e-8)return;
    const z=u*a[2]+v*b[2]+w*c[2];
    if(z>=-25&&z<=25&&z>nearest){nearest=z;hit={triangle,point:ids[[u,v,w].indexOf(Math.max(u,v,w))]};}
  });
  return hit;
}
const vertex=`
attribute vec3 position;
attribute float value;
attribute vec3 barycentric;
attribute float selected;
uniform vec2 angles;
uniform float zoom;
uniform float aspect;
uniform vec2 pan;
varying float v_value;
varying vec3 v_bary;
varying float v_selected;
void main(){
 float cy=cos(angles.x),sy=sin(angles.x),cx=cos(angles.y),sx=sin(angles.y);
 vec3 p=vec3(position.x*cy+position.z*sy,position.y,-position.x*sy+position.z*cy);
 p=vec3(p.x,p.y*cx-p.z*sx,p.y*sx+p.z*cx);
 gl_Position=vec4(p.x*zoom/aspect+pan.x,p.y*zoom+pan.y,-p.z*.04,1.0);
 v_value=value;v_bary=barycentric;v_selected=selected;
}`;
const fragment=`
precision mediump float;
varying float v_value;
varying vec3 v_bary;
varying float v_selected;
uniform bool wire;
uniform bool hasField;
uniform bool light;
vec3 ramp(float t){
 vec3 a=vec3(.15,.25,.47),b=vec3(.08,.57,.63),c=vec3(.55,.77,.60),d=vec3(.99,.75,.40),e=vec3(.85,.30,.18);
 float x=clamp(t,0.0,1.0)*4.0;
 if(x<1.0)return mix(a,b,x);if(x<2.0)return mix(b,c,x-1.0);if(x<3.0)return mix(c,d,x-2.0);return mix(d,e,x-3.0);
}
void main(){
 vec3 color=hasField?ramp(v_value):vec3(.33,.52,.62);
 float edge=min(v_bary.x,min(v_bary.y,v_bary.z));
 if(wire && edge<.012)color=mix(color,light?vec3(.10,.17,.23):vec3(.06,.12,.18),.7);
 if(v_selected>.5 && edge<.035)color=vec3(1.0,.93,.68);
 gl_FragColor=vec4(color,1.0);
}`;
export class MeshView {
  constructor(canvas,{onPick=()=>{},onError=()=>{}}={}) {
    this.canvas=canvas;this.onPick=onPick;this.onError=onError;this.options={};this.visible=false;this.ranges=new WeakMap();this.reset();
    this.gl=canvas.getContext('webgl',{alpha:true,depth:true,antialias:true});
    if(this.gl){try{this.initialize();}catch(e){this.gl=null;onError(e.message);}}
    if(!this.gl)this.software=canvas.getContext('2d');
    canvas.addEventListener('webglcontextlost',e=>{e.preventDefault();this.lost=true;onError('Graphics context lost. Data is retained; the viewport will recover when graphics returns.');});
    canvas.addEventListener('webglcontextrestored',()=>{try{this.lost=false;this.initialize();this.key=null;this.upload();this.invalidate();}catch(e){this.lost=true;onError(e.message);}});
    let drag=null;
    canvas.addEventListener('pointerdown',e=>{drag={id:e.pointerId,x:e.clientX,y:e.clientY,distance:0};canvas.setPointerCapture(e.pointerId);});
    canvas.addEventListener('pointermove',e=>{
      if(!drag||e.pointerId!==drag.id)return;
      const dx=e.clientX-drag.x,dy=e.clientY-drag.y;drag.distance+=Math.abs(dx)+Math.abs(dy);
      if(e.shiftKey){this.pan[0]+=dx/canvas.clientWidth*2;this.pan[1]-=dy/canvas.clientHeight*2;}
      else{this.yaw+=dx*.008;this.pitch=clamp(this.pitch+dy*.008,-1.55,1.55);}
      drag.x=e.clientX;drag.y=e.clientY;this.invalidate();
    });
    canvas.addEventListener('pointerup',e=>{if(!drag||e.pointerId!==drag.id)return;if(drag.distance<5)this.pick(e);drag=null;});
    canvas.addEventListener('pointercancel',()=>{drag=null;});
    canvas.addEventListener('wheel',e=>{e.preventDefault();this.zoomBy(Math.exp(-e.deltaY*.001));},{passive:false});
    canvas.addEventListener('keydown',e=>{
      if(!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','+','-','=','0'].includes(e.key))return;e.preventDefault();
      if(e.key==='0')this.reset();else if(['+','=','-'].includes(e.key))this.zoomBy(e.key==='-'?.9:1.1);
      else{this.yaw+=e.key==='ArrowLeft'?-.08:e.key==='ArrowRight'?.08:0;this.pitch=clamp(this.pitch+(e.key==='ArrowUp'?-.08:e.key==='ArrowDown'?.08:0),-1.55,1.55);this.invalidate();}
    });
    new ResizeObserver(()=>this.invalidate()).observe(canvas);
  }
  initialize() {
    const gl=this.gl,compile=(type,source)=>{
      const shader=gl.createShader(type);gl.shaderSource(shader,source);gl.compileShader(shader);
      if(!gl.getShaderParameter(shader,gl.COMPILE_STATUS)){const info=gl.getShaderInfoLog(shader);gl.deleteShader(shader);throw new Error(`Mesh shader failed: ${info}`);}return shader;
    };
    const vs=compile(gl.VERTEX_SHADER,vertex),fs=compile(gl.FRAGMENT_SHADER,fragment),p=gl.createProgram();
    gl.attachShader(p,vs);gl.attachShader(p,fs);gl.linkProgram(p);gl.deleteShader(vs);gl.deleteShader(fs);
    if(!gl.getProgramParameter(p,gl.LINK_STATUS))throw new Error('Mesh shader could not link.');
    this.program=p;this.buffer=gl.createBuffer();
    this.uniforms=Object.fromEntries(['angles','zoom','aspect','pan','wire','hasField','light'].map(k=>[k,gl.getUniformLocation(p,k)]));
    gl.useProgram(p);gl.bindBuffer(gl.ARRAY_BUFFER,this.buffer);
    for(const [name,size,offset] of [['position',3,0],['value',1,3],['barycentric',3,4],['selected',1,7]]) {
      const location=gl.getAttribLocation(p,name);gl.enableVertexAttribArray(location);gl.vertexAttribPointer(location,size,gl.FLOAT,false,32,offset*4);
    }
    gl.enable(gl.DEPTH_TEST);gl.depthFunc(gl.LEQUAL);gl.disable(gl.CULL_FACE);
  }
  reset(){
    this.yaw=-.45;this.pitch=.48;this.zoom=1.2;this.pan=[0,0];
    if(this.bounds&&this.options.mesh){
      const {center,span}=this.bounds,camera={yaw:this.yaw,pitch:this.pitch,zoom:1,aspect:1},aspect=Math.max(1,this.canvas.clientWidth)/Math.max(1,this.canvas.clientHeight);let x=0,y=0;
      for(const p of this.options.mesh.positions_m){const q=projectPoint(p.map((v,k)=>(v-center[k])/span*1.6),camera);x=Math.max(x,Math.abs(q[0]));y=Math.max(y,Math.abs(q[1]));}
      this.zoom=clamp(Math.min(.78*aspect/Math.max(x,1e-9),.72/Math.max(y,1e-9)),.1,8);
    }
    this.invalidate();
  }
  zoomBy(f){this.zoom=clamp(this.zoom*f,.1,8);this.invalidate();}
  set(options) {
    const changed=this.options.mesh!==options.mesh;this.options=options;this.visible=!!options.mesh;
    if(changed){this.softwareKey=null;this.lastError=null;this.bounds=options.mesh?boundsOf(options.mesh.positions_m):null;this.key=null;this.reset();}
    if(this.visible){this.upload();this.invalidate();}
  }
  range() {
    const {mesh,fieldIndex,component='magnitude'}=this.options,field=mesh?.fields[fieldIndex];if(!field)return [0,1];
    let ranges=this.ranges.get(field);if(!ranges){ranges=new Map();this.ranges.set(field,ranges);}
    if(!ranges.has(component))ranges.set(component,fieldRange(field,component));return ranges.get(component);
  }
  upload() {
    const {mesh,frame=0,fieldIndex=0,component='magnitude',deform=0,entity=0,displacementIndex=null}=this.options;
    if(this.lost||!mesh)return;
    const key=[frame,fieldIndex,component,deform,entity,displacementIndex].join(':');if(key===this.key)return;
    const field=mesh.fields[fieldIndex],displacement=displacementIndex!=null?mesh.fields[displacementIndex]:null;
    const {center,span}=this.bounds;
    const points=mesh.positions_m.map((p,i)=>p.map((v,k)=>((v-center[k])/span+(displacement?displacement.values[frame][i*3+k]/span*deform:0))*1.6));
    if(points.some(p=>p.some(v=>!Number.isFinite(v)||Math.abs(v)>20)))throw new Error('Deformation exceeds the viewport range. Reduce its display multiplier; source data is unchanged.');
    const [min,max]=this.range(),data=new Float32Array(mesh.triangles.length*3*8);let at=0;
    mesh.triangles.forEach((triangle,j)=>triangle.forEach((point,k)=>{
      const source=field?.association==='cell'?mesh.triangle_cells[j]:point;
      const value=field?fieldValue(field,frame,source,component):0;
      data.set([...points[point],max===min?.5:(value-min)/(max-min),k===0?1:0,k===1?1:0,k===2?1:0,source===entity?1:0],at);at+=8;
    }));
    this.points=points;this.key=key;this.vertices=data;
    if(!this.gl)return;
    this.gl.bindBuffer(this.gl.ARRAY_BUFFER,this.buffer);this.gl.bufferData(this.gl.ARRAY_BUFFER,data,this.gl.DYNAMIC_DRAW);
  }
  camera(){return {yaw:this.yaw,pitch:this.pitch,zoom:this.zoom,pan:this.pan,aspect:Math.max(1,this.canvas.clientWidth)/Math.max(1,this.canvas.clientHeight)};}
  pick(event) {
    const mesh=this.options.mesh;if(!mesh||!this.points)return;
    const r=this.canvas.getBoundingClientRect(),camera=this.camera();
    const points=this.points.map(p=>projectPoint(p,camera));
    const hit=pickTriangle(points,mesh.triangles,(event.clientX-r.left)/r.width*2-1,1-(event.clientY-r.top)/r.height*2);
    if(hit)this.onPick({...hit,cell:mesh.triangle_cells[hit.triangle]});
  }
  invalidate(){if(this.pending||!this.visible||document.hidden)return;this.pending=requestAnimationFrame(()=>{this.pending=null;this.draw();});}
  drawSoftware(){
    if(!this.software||!this.vertices)return;
    const width=this.canvas.clientWidth,height=this.canvas.clientHeight,scale=Math.min(1,760/Math.max(1,width),500/Math.max(1,height));
    const w=Math.max(1,Math.round(width*scale)),h=Math.max(1,Math.round(height*scale));
    const camera=this.camera(),key=JSON.stringify([this.key,camera,w,h,this.options.wire,document.documentElement.dataset.theme]);
    if(key===this.softwareKey)return;
    try{
      const pixels=rasterize(this.vertices,p=>projectPoint(p,camera),w,h,{wire:this.options.wire!==false,hasField:!!this.options.mesh.fields[this.options.fieldIndex],light:document.documentElement.dataset.theme==='light'});
      this.canvas.width=w;this.canvas.height=h;this.software.putImageData(new ImageData(pixels,w,h),0,0);this.softwareKey=key;
    }catch(e){this.canvas.width=w;this.canvas.height=h;if(this.lastError!==e.message){this.lastError=e.message;this.onError(e.message);}}
  }
  draw() {
    const gl=this.gl;if(this.lost||!this.visible)return;
    if(!gl){this.drawSoftware();return;}
    const dpr=Math.min(window.devicePixelRatio||1,2),w=Math.max(1,Math.round(this.canvas.clientWidth*dpr)),h=Math.max(1,Math.round(this.canvas.clientHeight*dpr));
    if(this.canvas.width!==w||this.canvas.height!==h){this.canvas.width=w;this.canvas.height=h;}
    gl.viewport(0,0,w,h);gl.clearColor(0,0,0,0);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);gl.useProgram(this.program);
    const u=this.uniforms;gl.uniform2f(u.angles,this.yaw,this.pitch);gl.uniform1f(u.zoom,this.zoom);gl.uniform1f(u.aspect,w/h);gl.uniform2fv(u.pan,this.pan);
    gl.uniform1i(u.wire,this.options.wire!==false);gl.uniform1i(u.hasField,!!this.options.mesh.fields[this.options.fieldIndex]);
    gl.uniform1i(u.light,document.documentElement.dataset.theme==='light');gl.drawArrays(gl.TRIANGLES,0,this.options.mesh.triangles.length*3);
  }
}
