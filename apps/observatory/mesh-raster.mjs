// Bounded, depth-tested fallback for hosts without WebGL. Same normalized vertex
// buffer, ramp and camera projection as the GPU; never changes source fields.
const stops=[[.15,.25,.47],[.08,.57,.63],[.55,.77,.60],[.99,.75,.40],[.85,.30,.18]];
export function ramp(t){const x=Math.max(0,Math.min(1,t))*4,i=Math.min(3,Math.floor(x)),f=x-i;return stops[i].map((v,k)=>v+(stops[i+1][k]-v)*f);}
export const MAX_PIXEL_VISITS=24000000;
export function rasterize(vertices,project,width,height,{wire=true,hasField=true,light=false}={}){
  if(!Number.isSafeInteger(width)||!Number.isSafeInteger(height)||width<1||height<1||width*height>400000)throw new Error('Software viewport exceeds its pixel budget.');
  const projected=new Float64Array(vertices.length/8*3);let visits=0;
  for(let i=0,j=0;i<vertices.length;i+=8,j+=3){const p=project([vertices[i],vertices[i+1],vertices[i+2]]);projected[j]=(p[0]+1)*width/2;projected[j+1]=(1-p[1])*height/2;projected[j+2]=p[2];}
  const boxes=[];
  // Preflight all coverage before publishing an image. Do not drop triangles.
  for(let i=0;i<projected.length;i+=9){
    const x0=Math.max(0,Math.floor(Math.min(projected[i],projected[i+3],projected[i+6]))),x1=Math.min(width-1,Math.ceil(Math.max(projected[i],projected[i+3],projected[i+6])));
    const y0=Math.max(0,Math.floor(Math.min(projected[i+1],projected[i+4],projected[i+7]))),y1=Math.min(height-1,Math.ceil(Math.max(projected[i+1],projected[i+4],projected[i+7])));
    visits+=Math.max(0,x1-x0+1)*Math.max(0,y1-y0+1);
    if(visits>MAX_PIXEL_VISITS)throw new Error('Mesh overlap exceeds the software viewport budget. Zoom out or use a WebGL-capable browser. Data and probes are retained.');
    boxes.push([x0,x1,y0,y1]);
  }
  const rgba=new Uint8ClampedArray(width*height*4),depth=new Float64Array(width*height).fill(-Infinity);
  for(let tri=0;tri<boxes.length;tri++){
    const p=tri*9,v=tri*24,[x0,x1,y0,y1]=boxes[tri];
    const ax=projected[p],ay=projected[p+1],bx=projected[p+3],by=projected[p+4],cx=projected[p+6],cy=projected[p+7];
    const det=(by-cy)*(ax-cx)+(cx-bx)*(ay-cy);if(Math.abs(det)<1e-12)continue;
    for(let y=y0;y<=y1;y++)for(let x=x0;x<=x1;x++){
      const u=((by-cy)*(x+.5-cx)+(cx-bx)*(y+.5-cy))/det;
      const b=((cy-ay)*(x+.5-cx)+(ax-cx)*(y+.5-cy))/det,c=1-u-b;
      if(Math.min(u,b,c)<0)continue;
      const z=u*projected[p+2]+b*projected[p+5]+c*projected[p+8],pixel=y*width+x;
      if(z<depth[pixel]||z>25||z< -25)continue;depth[pixel]=z;
      let color=hasField?ramp(u*vertices[v+3]+b*vertices[v+11]+c*vertices[v+19]):[.33,.52,.62];
      const edge=Math.min(u,b,c),selected=u*vertices[v+7]+b*vertices[v+15]+c*vertices[v+23];
      if(wire&&edge<.012)color=color.map((a,k)=>a*.3+(light?[.10,.17,.23]:[.06,.12,.18])[k]*.7);
      if(selected>.5&&edge<.035)color=[1,.93,.68];
      const at=pixel*4;rgba[at]=color[0]*255;rgba[at+1]=color[1]*255;rgba[at+2]=color[2]*255;rgba[at+3]=255;
    }
  }
  return rgba;
}
