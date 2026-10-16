// Mesh mode inside Observatory's existing run library, chart and export workflow.
import {MeshView} from './mesh-view.mjs';
import {fieldValue,frameAt,probe,boundsOf} from './mesh-data.mjs';
import {format,drawPlot,palette} from './visuals.mjs';
const $=id=>document.getElementById(id);
const node=(tag,text)=>{const e=document.createElement(tag);e.textContent=text;return e;};
export class MeshWorkspace {
  constructor({notify,onBack}) {
    this.notify=notify;this.active=false;this.frame=0;this.field=0;this.entity=0;this.component='magnitude';
    const inspector=document.createElement('div');inspector.id='mesh-inspector';
    // Static application markup only; all file-supplied labels use textContent.
    inspector.innerHTML=`<div class="panel-heading"><span class="eyebrow">01 / FIELD EXPLORER</span><span class="tag">RECORDED DATA</span></div>
      <h2>Inspect the continuum</h2><p class="small muted">One mesh. Every physical field.</p>
      <div class="mesh-source-card"><span class="eyebrow">SOURCE MESH</span><strong id="mesh-size"></strong><p class="small muted">Surface geometry is displayed to scale. No solver runs in this view.</p></div>
      <label for="mesh-field">Physical field</label><select id="mesh-field"></select>
      <div id="mesh-component-box"><label for="mesh-component">Vector display</label><select id="mesh-component"><option value="magnitude">Magnitude |v|</option><option value="0">X component</option><option value="1">Y component</option><option value="2">Z component</option></select></div>
      <label for="mesh-entity" id="mesh-entity-label">Probe entity (zero-based)</label><input id="mesh-entity" type="number" min="0" step="1" value="0">
      <p class="small muted mesh-hint">Click a visible element or enter an original entity index. Values below are source samples, not color estimates.</p>
      <label class="mesh-check"><input id="mesh-wire" type="checkbox"> Show element edges</label>
      <div id="mesh-deformation-box"><label for="mesh-displacement-field">Displacement vector</label><select id="mesh-displacement-field"><option value="">Reference geometry only</option></select><label for="mesh-deformation">Displacement display multiplier <span>×</span></label><input id="mesh-deformation" type="number" min="0" max="100" step="any" value="0"><p class="small muted">0 = reference geometry. Scaling is a visual aid and never alters stored values.</p></div>
      <div class="mesh-source-card"><span class="eyebrow">TRUST THE SCALE</span><p class="small muted">Color limits span every recorded frame. Cell values stay discontinuous across material boundaries. Vertex colors use linear interpolation for display only.</p></div>
      <button id="mesh-new-experiment" class="button secondary">← New native experiment</button>`;
    document.querySelector('.inspector').append(inspector);
    const canvas=document.createElement('canvas');canvas.id='mesh-scene';canvas.tabIndex=0;
    canvas.setAttribute('aria-label','Recorded surface mesh. Drag to rotate, Shift-drag to pan, scroll to zoom. Click to probe; arrow keys rotate.');
    document.querySelector('.canvas-stage').prepend(canvas);
    const footer=document.createElement('div');footer.id='mesh-legend';
    footer.innerHTML='<div><strong id="mesh-field-caption"></strong><span id="mesh-association"></span></div><div class="mesh-scale"><div class="mesh-ramp"></div><div class="mesh-scale-labels"><span id="mesh-min"></span><span id="mesh-max"></span></div></div>';
    document.querySelector('.scene-footer').append(footer);
    const timeline=document.createElement('div');timeline.id='mesh-timeline';
    timeline.innerHTML=`<div class="mesh-transport"><button id="mesh-prev" class="small-button" aria-label="Previous recorded frame">|‹</button><button id="mesh-play" class="small-button" aria-label="Play recorded frames">▶ Play</button><button id="mesh-next" class="small-button" aria-label="Next recorded frame">›|</button><output id="mesh-clock"></output><label class="sr-only" for="mesh-speed">Playback speed</label><select id="mesh-speed"><option value="0.1">0.1×</option><option value="0.25">0.25×</option><option value="1" selected>1×</option><option value="2">2×</option></select></div><label class="sr-only" for="mesh-frame">Recorded frame</label><input id="mesh-frame" type="range" min="0" max="0" value="0"><span id="mesh-frame-note" class="small muted"></span>`;
    document.querySelector('.scene-panel').append(timeline);
    const link=document.createElement('link');link.rel='stylesheet';link.href='/mesh.css';document.head.append(link);
    this.originalEmpty=$('empty-chart').innerHTML;
    this.originalText=new Map();
    for(const selector of ['.metric>span','.metric>small','#signal-heading','.selection-card>.eyebrow']) for(const el of document.querySelectorAll(selector))this.originalText.set(el,el.textContent);
    this.view=new MeshView(canvas,{onPick:hit=>{const f=this.mesh.fields[this.field];this.entity=f?.association==='cell'?hit.cell:hit.point;$('mesh-entity').value=this.entity;this.render();},onError:m=>notify(m,true)});
    $('mesh-field').addEventListener('change',()=>{this.field=Number($('mesh-field').value);this.entity=0;this.component='magnitude';$('mesh-component').value='magnitude';this.render();});
    $('mesh-component').addEventListener('change',()=>{this.component=$('mesh-component').value;this.render();});
    $('mesh-entity').addEventListener('change',()=>{
      const id=Number($('mesh-entity').value),count=this.mesh.fields[this.field]?.association==='cell'?this.mesh.cell_count:this.mesh.positions_m.length;
      if(!Number.isSafeInteger(id)||id<0||id>=count){notify(`Entity must be an integer from 0 to ${count-1}.`,true);$('mesh-entity').value=this.entity;return;}
      this.entity=id;this.render();
    });
    $('mesh-displacement-field').addEventListener('change',()=>{this.displacement=$('mesh-displacement-field').value===''?null:Number($('mesh-displacement-field').value);this.deform=0;$('mesh-deformation').value='0';this.render();});
    $('mesh-deformation').addEventListener('change',()=>{
      const value=Number($('mesh-deformation').value);if(!Number.isFinite(value)||value<0||value>100){$('mesh-deformation').value=this.deform;notify('Display multiplier must be between 0 and 100.',true);return;}
      const old=this.deform;this.deform=value;
      try{this.updateGPU();}catch(e){this.deform=old;$('mesh-deformation').value=old;notify(e.message,true);}this.render();
    });
    $('mesh-wire').addEventListener('change',()=>this.render());
    $('mesh-frame').addEventListener('input',()=>{this.pause();this.setFrame(Number($('mesh-frame').value));});
    $('mesh-prev').addEventListener('click',()=>{this.pause();this.setFrame(this.frame-1);});
    $('mesh-next').addEventListener('click',()=>{this.pause();this.setFrame(this.frame+1);});
    $('mesh-play').addEventListener('click',()=>this.playing?this.pause():this.play());
    $('mesh-speed').addEventListener('change',()=>{if(this.playing){this.pause();this.play();}});
    $('mesh-new-experiment').addEventListener('click',()=>{this.pause();onBack();});
    document.addEventListener('visibilitychange',()=>{if(document.hidden)this.pause();else this.view.invalidate();});
    window.addEventListener('pagehide',()=>this.pause());
    this.setVisible(false);
  }
  setVisible(active){this.active=active;document.body.classList.toggle('mesh-mode',active);for(const id of ['mesh-inspector','mesh-scene','mesh-legend','mesh-timeline'])$(id).hidden=!active;}
  enter(run,baseline) {
    this.setVisible(true);this.run=run;this.baseline=baseline;
    if(this.mesh!==run.snapshot.mesh){
      this.pause();this.mesh=run.snapshot.mesh;this.frame=0;this.field=0;this.entity=0;this.component='magnitude';this.deform=0;this.comparison=null;
      this.displacement=this.mesh.displacement_field??null;
      $('mesh-displacement-field').replaceChildren();const reference=node('option','Reference geometry only');reference.value='';$('mesh-displacement-field').append(reference);
      this.mesh.fields.forEach((f,i)=>{if(f.association==='point'&&f.components===3&&f.unit==='m'){const option=node('option',f.label+' · m');option.value=i;$('mesh-displacement-field').append(option);}});
      $('mesh-displacement-field').value=this.displacement??'';
      $('mesh-deformation').value='0';$('mesh-component').value='magnitude';$('mesh-field').replaceChildren();
      this.mesh.fields.forEach((f,i)=>{const option=node('option',`${f.label} · ${f.unit} · ${f.association}`);option.value=i;$('mesh-field').append(option);});
      if(!this.mesh.fields.length)$('mesh-field').append(node('option','Geometry only · no fields'));
    }
    this.render();
  }
  exit(){if(!this.active)return;this.pause();this.setVisible(false);this.view.set({mesh:null});for(const [el,text] of this.originalText)el.textContent=text;$('empty-chart').innerHTML=this.originalEmpty;}
  updateGPU(){this.view.set({mesh:this.mesh,frame:this.frame,fieldIndex:this.field,component:this.component,deform:this.deform,displacementIndex:this.displacement,entity:this.entity,wire:$('mesh-wire').checked});}
  trace(){return this.mesh?probe(this.mesh,this.field,this.entity,this.component):null;}
  baselineTrace(){
    const other=this.baseline?.snapshot?.mesh,field=this.mesh.fields[this.field];
    if(!other||other===this.mesh||!field)return null;
    const key=[this.field,this.entity,this.component].join(':');
    if(this.comparison?.other===other&&this.comparison.key===key)return this.comparison.trace;
    const topology=m=>JSON.stringify([m.positions_m,m.triangles,m.triangle_cells,m.cell_count,m.source_cells,m.source_types]);
    const index=other.fields.findIndex(f=>f.label===field.label&&f.unit===field.unit&&f.components===field.components&&f.association===field.association);
    const identityKnown=field.association==='point'||!!(this.mesh.source_cells&&other.source_cells);
    const trace=index>=0&&identityKnown&&topology(other)===topology(this.mesh)?probe(other,index,this.entity,this.component):null;
    this.comparison={other,key,trace};return trace;
  }
  render(){
    if(!this.active)return;
    const m=this.mesh,f=m.fields[this.field],t=m.time_s[this.frame],bounds=boundsOf(m.positions_m);
    try{this.updateGPU();}catch(e){
      // A valid recorded field may exceed an exaggerated display scale in a later frame.
      // Fall back to undeformed geometry; never lose the selected field or its values.
      this.pause();this.deform=0;$('mesh-deformation').value='0';this.updateGPU();this.notify(e.message,true);
    }
    $('scene-heading').textContent='Recorded field explorer';$('scene-state').textContent=`${this.run.origin.toUpperCase()} · ${m.time_s.length>1?'RECORDED SEQUENCE':'STATIC SNAPSHOT'}`;
    $('scene-caption').textContent=this.view.gl?'GEOMETRY TO SCALE · ORIGINAL SOURCE ENTITIES':'SOFTWARE VIEWPORT · GEOMETRY TO SCALE · SOURCE PROBES';
    $('scene-subtitle').textContent=`Drag to orbit · Shift-drag to pan · click to probe${this.deform ? ` · displacement ×${this.deform} (visual only)` : ' · reference geometry'}`;
    $('mesh-size').textContent=`${m.positions_m.length.toLocaleString()} points / ${m.cell_count.toLocaleString()} cells`;
    $('mesh-field').disabled=!f;$('mesh-component-box').hidden=f?.components!==3;$('mesh-entity').disabled=!f;
    $('mesh-component').value=this.component;$('mesh-field').value=String(this.field);
    $('mesh-entity').value=this.entity;$('mesh-entity').max=(f?.association==='cell'?m.cell_count:m.positions_m.length)-1;
    $('mesh-entity-label').textContent=`Probe ${f?.association||'point'} index (zero-based)`;
    $('mesh-deformation-box').hidden=$('mesh-displacement-field').options.length<2;$('mesh-deformation').disabled=this.displacement==null;
    const range=this.view.range();$('mesh-min').textContent=f?`${format(range[0],5)} ${f.unit}`:'—';$('mesh-max').textContent=f?`${format(range[1],5)} ${f.unit}`:'—';
    $('mesh-field-caption').textContent=f?f.label:'Geometry only';$('mesh-association').textContent=f?`${f.association} data · all-frame scale${f.components===3?' · '+(this.component==='magnitude'?'magnitude':'XYZ'[Number(this.component)]):''}`:'No physical values inferred';
    $('mesh-frame').max=m.time_s.length-1;$('mesh-frame').value=this.frame;$('mesh-frame').disabled=m.time_s.length===1;
    $('mesh-clock').textContent=m.time_s.length===1?'Static snapshot':`${format(t,6)} s`;
    $('mesh-frame-note').textContent=`Frame ${this.frame+1} / ${m.time_s.length} · recorded samples only · no time interpolation`;
    $('mesh-play').disabled=m.time_s.length<2;$('mesh-prev').disabled=this.frame===0;$('mesh-next').disabled=this.frame===m.time_s.length-1;
    let min=Infinity,max=-Infinity;const n=f?f.values[this.frame].length/f.components:0;
    for(let i=0;i<n;i++){const v=fieldValue(f,this.frame,i,this.component);min=Math.min(min,v);max=Math.max(max,v);}
    const metricValues=[f?min:null,f?max:null,f?fieldValue(f,this.frame,this.entity,this.component):null,m.time_s.length>1?t:null];
    const names=['Field minimum','Field maximum','Probe value','Recorded time'];
    const notes=['All source entities · this frame','All source entities · this frame',`${f?.association||'Entity'} ${this.entity} · exact source sample`,m.time_s.length>1?'Playback is a recorded-frame hold':'No solver time asserted'];
    [...document.querySelectorAll('.metric')].forEach((card,i)=>{card.querySelector('span').textContent=names[i];const strong=card.querySelector('strong');strong.replaceChildren(node('span',format(metricValues[i],5)));if(metricValues[i]!=null)strong.append(node('em',i===3?'s':f.unit));card.querySelector('small').textContent=notes[i];});
    document.querySelector('.selection-card>.eyebrow').textContent='EXACT SOURCE PROBE';$('selected-layer').textContent=`${f?.association==='cell'?'Cell':'Point'} ${this.entity}`;
    const dl=$('layer-details');dl.replaceChildren();
    const rows=[['Association',f?.association||'No field'],['Field',f?.label||'Not supplied'],['Value',f?`${format(metricValues[2],7)} ${f.unit}`:'—'],['Displayed triangles',String(m.triangles.length)],['Reference extent (m)',bounds.max.map((v,i)=>format(v-bounds.min[i],5)).join(' × ')]];
    if(f?.association==='point')rows.splice(3,0,['Reference XYZ (m)',m.positions_m[this.entity].map(v=>format(v,7)).join(', ')]);
    if(f?.components===3)rows.splice(3,0,['Original XYZ components',f.values[this.frame].slice(this.entity*3,this.entity*3+3).map(v=>format(v,7)).join(', ')]);
    rows.push(['Source',this.run.snapshot.source.material]);for(const [k,v] of rows)dl.append(node('dt',k),node('dd',v));
    $('signal-heading').textContent='Probe history · recorded samples';$('export-csv').disabled=!f;$('export-json').disabled=false;
    $('listen-cold').disabled=true;$('listen-hot').disabled=true;this.drawProbe();
  }
  drawProbe(){
    if(!this.active)return;
    const trace=this.trace(),m=this.mesh,colors=palette(),lines=[];let unit='';
    if(trace&&m.time_s.length>1){const s=trace.series[0];unit=s.unit;lines.push({x:trace.time_s,y:s.values,color:colors[0]});}
    const comparison=this.baselineTrace();if(comparison&&comparison.time_s.length>1&&m.time_s.length>1)lines.push({x:comparison.time_s,y:comparison.series[0].values,color:colors[2],dashed:true});
    const range=m.time_s.length>1?[m.time_s[0],m.time_s.at(-1)]:[0,1];
    this.plot=drawPlot($('trace'),lines,{range,cursor:m.time_s.length>1?m.time_s[this.frame]:null,xunit:'s',yunit:unit});
    $('empty-chart').hidden=!!lines.length;
    $('empty-chart').textContent=!trace?'No field values supplied.':'One recorded snapshot. A time history is not inferred.';
    $('channels').replaceChildren();if(trace)$('channels').append(node('span',`${trace.series[0].label}: ${format(trace.series[0].values[this.frame],7)} ${trace.series[0].unit}`));
    if(this.baseline&&this.baseline!==this.run)$('channels').append(node('span',comparison?'┄ Pinned run · same source geometry and field':'Pinned run is incompatible with this field/entity geometry.'));
  }
  setFrame(index){this.frame=Math.max(0,Math.min(this.mesh.time_s.length-1,index));this.render();}
  cursor(time){this.pause();this.setFrame(frameAt(this.mesh.time_s,time));}
  step(direction){this.pause();this.setFrame(this.frame+direction);}
  play(){
    if(!this.active||this.mesh.time_s.length<2)return;
    if(this.frame===this.mesh.time_s.length-1)this.setFrame(0);
    const start=performance.now(),time=this.mesh.time_s[this.frame],speed=Number($('mesh-speed').value);
    this.playing=true;$('mesh-play').textContent='Ⅱ Pause';$('mesh-play').setAttribute('aria-label','Pause recorded frames');
    const tick=now=>{
      if(!this.playing)return;const target=time+(now-start)/1000*speed,index=frameAt(this.mesh.time_s,target);
      if(index!==this.frame)this.setFrame(index);
      if(target>=this.mesh.time_s.at(-1)){this.pause();return;}this.timer=requestAnimationFrame(tick);
    };this.timer=requestAnimationFrame(tick);
  }
  pause(){this.playing=false;if(this.timer)cancelAnimationFrame(this.timer);this.timer=null;if($('mesh-play')){$('mesh-play').textContent='▶ Play';$('mesh-play').setAttribute('aria-label','Play recorded frames');}}
  screenshot(download,name){
    this.view.draw();const source=this.view.canvas,canvas=document.createElement('canvas');canvas.width=source.width;canvas.height=source.height+100;
    const c=canvas.getContext('2d'),css=getComputedStyle(document.documentElement);c.fillStyle=css.getPropertyValue('--panel');c.fillRect(0,0,canvas.width,canvas.height);c.drawImage(source,0,0);
    c.fillStyle=css.getPropertyValue('--text');c.font='14px monospace';const f=this.mesh.fields[this.field],range=this.view.range();
    c.fillText(`${this.run.origin} · ${this.run.snapshot.title}`.slice(0,100),16,source.height+24);
    c.fillText(f?`${f.label} (${f.unit}) · ${f.association} · full-history range ${format(range[0],5)} to ${format(range[1],5)}`:'Geometry only',16,source.height+49);
    c.fillText(`Frame ${this.frame+1}/${this.mesh.time_s.length} · ${this.mesh.time_s.length>1?this.mesh.time_s[this.frame]+' s':'static'} · displacement ×${this.deform} (display only)`,16,source.height+74);
    canvas.toBlob(blob=>{if(blob)download(blob,'image/png',name);});
  }
}
