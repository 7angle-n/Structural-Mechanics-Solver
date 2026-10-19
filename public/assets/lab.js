/* Structural Mechanics Solver homepage: homepage visuals. The hero/CTA scene is a REAL 3D truss solved here by the direct stiffness method. */
const q=s=>document.querySelector(s),qa=s=>[...document.querySelectorAll(s)],f2=(v,n=2)=>(+v.toFixed(n)).toString();
const RM=matchMedia("(prefers-reduced-motion:reduce)").matches;
function gauss(A,b){const n=b.length;A=A.map((r,i)=>[...r,b[i]]);for(let c=0;c<n;c++){let p=c;for(let r=c+1;r<n;r++)if(Math.abs(A[r][c])>Math.abs(A[p][c]))p=r;if(Math.abs(A[p][c])<1e-9)throw Error("singular");[A[c],A[p]]=[A[p],A[c]];for(let r=0;r<n;r++)if(r!=c){const k=A[r][c]/A[c][c];for(let j=c;j<=n;j++)A[r][j]-=k*A[c][j]}}return A.map((r,i)=>r[n]/r[i])}
/* ---- model: 2-bay braced tower, pinned base, kN & m ---- */
const EA=200e6*0.002,ND=[],MB=[];
for(let z=0;z<2;z++)for(let iy=0;iy<2;iy++)for(let ix=0;ix<3;ix++)ND.push([ix*4,iy*3,z*4]);
const nid=(ix,iy,z)=>z*6+iy*3+ix;
for(let k=0;k<6;k++)MB.push([k,k+6]);
for(let iy=0;iy<2;iy++)for(let ix=0;ix<2;ix++)MB.push([nid(ix,iy,1),nid(ix+1,iy,1)]);
for(let ix=0;ix<3;ix++)MB.push([nid(ix,0,1),nid(ix,1,1)]);
for(let iy=0;iy<2;iy++)MB.push([nid(0,iy,0),nid(1,iy,1)],[nid(2,iy,0),nid(1,iy,1)]);
MB.push([nid(0,0,0),nid(0,1,1)],[nid(2,0,0),nid(2,1,1)],[nid(0,0,1),nid(1,1,1)],[nid(2,0,1),nid(1,1,1)]);
const LOADS=[[9,[40,0,0]],[7,[0,0,-30]],[10,[0,0,-30]]];
function solveTruss(load){const n=ND.length,K=Array.from({length:3*n},()=>Array(3*n).fill(0));
 MB.forEach(([i,j])=>{const d=[0,1,2].map(a=>ND[j][a]-ND[i][a]),L=Math.hypot(...d),c=d.map(v=>v/L),k=EA/L;for(let a=0;a<3;a++)for(let b=0;b<3;b++){const v=k*c[a]*c[b];K[3*i+a][3*i+b]+=v;K[3*j+a][3*j+b]+=v;K[3*i+a][3*j+b]-=v;K[3*j+a][3*i+b]-=v}});
 const fr=[];for(let i=6;i<12;i++)for(let a=0;a<3;a++)fr.push(3*i+a);
 const F=Array(3*n).fill(0);load.forEach(([i,f])=>f.forEach((v,a)=>F[3*i+a]+=v));
 const u0=gauss(fr.map(r=>fr.map(c=>K[r][c])),fr.map(r=>F[r])),u=Array(3*n).fill(0);fr.forEach((r,k)=>u[r]=u0[k]);
 const Nf=MB.map(([i,j])=>{const d=[0,1,2].map(a=>ND[j][a]-ND[i][a]),L=Math.hypot(...d);return EA/L*d.reduce((s,v,a)=>s+v/L*(u[3*j+a]-u[3*i+a]),0)});
 const um=Math.max(...[6,7,8,9,10,11].map(i=>Math.hypot(u[3*i],u[3*i+1],u[3*i+2])));return{u,Nf,um,Nm:Math.max(...Nf.map(Math.abs))}}
let SOL;try{SOL=solveTruss(LOADS)}catch(e){SOL=null}
const sm=x=>x<=0?0:x>=1?1:x*x*(3-2*x),HX=h=>[1,3,5].map(i=>parseInt(h.substr(i,2),16)),mixc=(a,b,t)=>`rgb(${HX(a).map((v,i)=>Math.round(v+(HX(b)[i]-v)*t)).join(",")})`;
/* ---- scene ---- */
function Scene(cv,o={}){
 const ctx=cv.getContext("2d"),S={yaw:-.55,pitch:.42,zoom:o.zoom||1,L:{loads:1,supports:1,undef:1,def:1,force:1,labels:0,...o.layers},drag:0,run:1,t:RM||o.still?9:0,c:0,user:0};
 const mag=SOL&&SOL.um?.55/SOL.um:0;
 function frame(){
  const w=cv.clientWidth,h=cv.clientHeight,d=devicePixelRatio||1;if(!w)return;if(cv.width!=Math.round(w*d)){cv.width=Math.round(w*d);cv.height=Math.round(h*d)}
  ctx.setTransform(d,0,0,d,0,0);ctx.clearRect(0,0,w,h);
  const T=S.t,k=Math.min(w*(o.wf||.06),h*.1)*S.zoom,ox=w*(o.cx||.5),oy=h*(o.cy||.56),cy=Math.cos(S.yaw),sy=Math.sin(S.yaw),cp=Math.cos(S.pitch),sp=Math.sin(S.pitch);
  const fade=o.still?1:T<.3?T/.3:T>15.4?Math.max(0,(16-T)/.6):1,lam=o.still||RM?1:T<13?sm((T-4.6)/3.5):1-sm((T-13)/2.5),la=o.still?1:sm((T-3.4)/1),na=o.still?1:sm((T-2.4)/.8);
  const rot=p=>{const x=p[0]-4,y=p[1]-1.5,z=p[2]-2,xr=x*cy-y*sy,yr=x*sy+y*cy;return[xr,z*cp+yr*sp,yr*cp-z*sp]};
  const pr=p=>{const r=rot(p),f=26/(26+r[2]);return[ox+r[0]*k*f,oy-r[1]*k*f]};
  const dp=i=>SOL?ND[i].map((v,a)=>v+SOL.u[3*i+a]*mag*lam):ND[i];
  ctx.globalAlpha=fade;ctx.lineCap="round";ctx.font="500 11px 'CM Typewriter','Latin Modern',monospace";
  ctx.strokeStyle="#12223f";ctx.lineWidth=1;for(let i=-1;i<=9;i+=2){let a=pr([i,-1.5,0]),b=pr([i,4.5,0]);ctx.beginPath();ctx.moveTo(...a);ctx.lineTo(...b);ctx.stroke()}for(let j=-1;j<=4;j+=1.5){let a=pr([-1,j,0]),b=pr([9,j,0]);ctx.beginPath();ctx.moveTo(...a);ctx.lineTo(...b);ctx.stroke()}
  const seg=(a,b,t)=>{const A=pr(a),B=pr(b);ctx.beginPath();ctx.moveTo(...A);ctx.lineTo(A[0]+(B[0]-A[0])*t,A[1]+(B[1]-A[1])*t);ctx.stroke()};
  MB.forEach(([i,j],m)=>{const p=o.still?1:sm((T-.1*m)/.7);if(p<=0)return;
   if(S.L.undef||!S.L.def||lam<.02){ctx.strokeStyle=S.L.def&&lam>.02?"#3b5f8a":"#7f95b8";ctx.lineWidth=S.L.def&&lam>.02?1:2;ctx.setLineDash(S.L.def&&lam>.02?[4,4]:[]);seg(ND[i],ND[j],p);ctx.setLineDash([])}
   if(S.L.def&&lam>.02&&p>=1){const f=SOL?SOL.Nf[m]:0;ctx.strokeStyle=S.L.force&&SOL?mixc("#7f95b8",f>0?"#38bdf8":"#fb7185",Math.min(1,Math.abs(f)/SOL.Nm*lam*1.2)):"#38bdf8";ctx.lineWidth=2.6;seg(dp(i),dp(j),1)}});
  if(S.L.supports)for(let i=0;i<6;i++){const s=pr(ND[i]);ctx.globalAlpha=fade*na;ctx.strokeStyle="#c9d6e8";ctx.lineWidth=1.4;ctx.beginPath();ctx.moveTo(s[0],s[1]);ctx.lineTo(s[0]-8,s[1]+13);ctx.lineTo(s[0]+8,s[1]+13);ctx.closePath();ctx.moveTo(s[0]-11,s[1]+16);ctx.lineTo(s[0]+11,s[1]+16);ctx.stroke()}
  ctx.globalAlpha=fade*na;ctx.fillStyle="#e6edf7";for(let i=0;i<12;i++){const s=pr(S.L.def&&lam>.02?dp(i):ND[i]);ctx.beginPath();ctx.arc(s[0],s[1],2.6,0,7);ctx.fill();if(S.L.labels){ctx.fillStyle="#8da2bf";ctx.fillText("N"+(i+1),s[0]+5,s[1]-5);ctx.fillStyle="#e6edf7"}}
  if(S.L.loads&&la>0){ctx.globalAlpha=fade;LOADS.forEach(([i,fv])=>{const mg=Math.hypot(...fv),dr=fv.map(v=>v/mg),H=ND[i],Tl=H.map((v,a)=>v-dr[a]*2.6*la),A=pr(Tl),B=pr(S.L.def&&lam>.02?dp(i):H);ctx.strokeStyle=ctx.fillStyle="#fbbf24";ctx.lineWidth=2.4;ctx.beginPath();ctx.moveTo(...A);ctx.lineTo(...B);ctx.stroke();const an=Math.atan2(B[1]-A[1],B[0]-A[0]);ctx.beginPath();ctx.moveTo(...B);ctx.lineTo(B[0]-10*Math.cos(an-.4),B[1]-10*Math.sin(an-.4));ctx.lineTo(B[0]-10*Math.cos(an+.4),B[1]-10*Math.sin(an+.4));ctx.fill();if(la>.6)ctx.fillText(mg+" kN",A[0]+6,A[1]-6)})}
  ctx.globalAlpha=fade*na;ctx.strokeStyle=ctx.fillStyle="#8da2bf";ctx.lineWidth=1;
  if(!o.compact){const dim=(a,b,t)=>{const A=pr(a),B=pr(b);ctx.beginPath();ctx.moveTo(...A);ctx.lineTo(...B);ctx.stroke();ctx.fillText(t,(A[0]+B[0])/2+6,(A[1]+B[1])/2+14)};dim([0,-.9,0],[8,-.9,0],"8.0 m");dim([8.9,0,0],[8.9,0,4],"4.0 m")}
  ctx.globalAlpha=1;if(!o.compact){const bx=44,by=h-58;[[1,0,0,"X","#fb7185"],[0,1,0,"Y","#34d399"],[0,0,1,"Z","#38bdf8"]].forEach(([x,y,z,n,c])=>{const xr=x*cy-y*sy,yr=x*sy+y*cy;ctx.strokeStyle=ctx.fillStyle=c;ctx.beginPath();ctx.moveTo(bx,by);ctx.lineTo(bx+xr*26,by-(z*cp+yr*sp)*26);ctx.stroke();ctx.fillText(n,bx+xr*32-3,by-(z*cp+yr*sp)*32+4)})}
  if(o.ro&&SOL){const f=lam;o.ro.textContent=`load factor ${f2(f)} · u_max ${f2(SOL.um*f*1000)} mm · N_max ${f2(SOL.Nm*f)} kN · deflections ×${Math.round(mag)}`}
 }
 let last=performance.now();function tick(ts){requestAnimationFrame(tick);if(!S.run||document.hidden)return;const dt=Math.min(.05,(ts-last)/1000);last=ts;
  if(!RM){S.c+=dt;if(!o.still)S.t=(S.t+dt)%16;if(!S.user)S.yaw=o.still?S.yaw+dt*.12:-.55+.45*Math.sin(S.c*.3)}frame()}
 requestAnimationFrame(tick);if("IntersectionObserver"in window)new IntersectionObserver(e=>S.run=e[0].isIntersecting).observe(cv);
 if(o.interactive){let lx,ly;cv.style.cursor="grab";cv.style.touchAction="pan-y";cv.addEventListener("pointerdown",e=>{S.drag=1;S.user=1;lx=e.clientX;ly=e.clientY;cv.setPointerCapture(e.pointerId)});
  cv.addEventListener("pointermove",e=>{if(!S.drag)return;S.yaw+=(e.clientX-lx)*.008;S.pitch=Math.max(.05,Math.min(1.35,S.pitch+(e.clientY-ly)*.006));lx=e.clientX;ly=e.clientY;if(RM)frame()});cv.addEventListener("pointerup",()=>S.drag=0)}
 frame();return{S,frame}}
/* ---- hero + CTA ---- */
const hero=Scene(q("#hero-cv"),{interactive:1,cx:innerWidth>820?.68:.5,wf:.05,ro:q("#hero-ro")});
qa("#hero-ctl [data-l]").forEach(b=>b.onclick=()=>{const l=b.dataset.l;hero.S.L[l]=hero.S.L[l]?0:1;b.classList.toggle("on",!!hero.S.L[l]);if(RM)hero.frame()});
qa("#hero-ctl [data-z]").forEach(b=>b.onclick=()=>{hero.S.zoom=Math.max(.6,Math.min(2,hero.S.zoom*(b.dataset.z=="+"?1.15:.87)));if(RM)hero.frame()});
Scene(q("#cta-cv"),{still:1,compact:1,layers:{labels:0},wf:.07,cx:.5});
/* ---- playground: closed-form beams + scaled truss solution ---- */
(function(){const $s=q("#pg-sys"),$p=q("#pg-p"),$a=q("#pg-a"),svg=q("#pg-svg"),out=q("#pg-out"),EI=40000,L=6,X=x=>40+x/L*520;
 const path=(xs,ys,y0,sc,col)=>`<path d="M${X(0)} ${y0} ${xs.map((x,i)=>`L${X(x)} ${y0-ys[i]*sc}`).join(" ")} L${X(L)} ${y0}Z" fill="${col}33" stroke="${col}" stroke-width="1.6"/>`;
 const tri=x=>`<path d="M${x} 84l-9 15h18z" fill="none" stroke="#c9d6e8" stroke-width="1.5"/>`,arr=(x,t,b)=>`<path d="M${x} ${t}V${b}m-5 -9l5 9l5 -9" fill="none" stroke="#fbbf24" stroke-width="2.4"/>`;
 function run(){const sys=$s.value,p=+$p.value,a=+$a.value;q("#pg-pv").textContent=p+" kN";q("#pg-av").textContent=a+" m";$a.disabled=sys!="ss";let g="",o="";
  if(sys=="truss"){const lm=p/40,mg=SOL.um?60/(SOL.um*Math.max(lm,.4)):0,P=(i,d)=>{const n=ND[i],u=d?[SOL.u[3*i]*lm,SOL.u[3*i+2]*lm]:[0,0];return[60+(n[0]+u[0]*(d?mg:0))/8*440,230-(n[2]+u[1]*(d?mg:0))/4*170]};
   MB.forEach(([i,j],m)=>{const A=P(i,0),B=P(j,0),C=P(i,1),D=P(j,1),f=SOL.Nf[m]*lm,t=Math.min(1,Math.abs(f)/(SOL.Nm*lm||1));g+=`<path d="M${A}L${B}" stroke="#3b5f8a" stroke-dasharray="3 3" fill="none"/><path d="M${C}L${D}" stroke="${mixc("#7f95b8",f>0?"#38bdf8":"#fb7185",t)}" stroke-width="2.4" fill="none"/>`});
   for(let i=0;i<3;i++)g+=`<path d="M${60+i*220} 232l-9 14h18z" fill="none" stroke="#c9d6e8"/>`;g+=`<path d="M16 40H56m-9 -5l9 5l-9 5" fill="none" stroke="#fbbf24" stroke-width="2.4"/><text x="14" y="30" fill="#fbbf24" font-size="12">${f2(40*lm,1)} kN</text>`;
   o=`<div>u_max<b>${f2(SOL.um*lm*1000)} mm</b></div><div>N_max<b>${f2(SOL.Nm*lm)} kN</b></div><div>Members<b>${MB.length}</b></div>`}
  else{const n=121,xs=[...Array(n).keys()].map(i=>i*L/(n-1));let V,M,Y,info;
   if(sys=="ss"){const b=L-a,RA=p*b/L,RB=p*a/L;V=xs.map(x=>x<a?RA:-RB);M=xs.map(x=>x<=a?RA*x:RB*(L-x));Y=xs.map(x=>x<=a?-p*b*x*(L*L-b*b-x*x)/(6*EI*L):-p*a*(L-x)*(L*L-a*a-(L-x)**2)/(6*EI*L));info=[["R_A",RA+" kN"],["R_B",RB+" kN"],["M_max",p*a*b/L+" kN·m"]]}
   else{V=xs.map(()=>p);M=xs.map(x=>-p*(L-x));Y=xs.map(x=>-p*x*x*(3*L-x)/(6*EI));info=[["R",p+" kN"],["M_wall",-p*L+" kN·m"]]}
   const dm=Math.max(...Y.map(Math.abs)),vm=Math.max(...V.map(Math.abs)),mm=Math.max(...M.map(Math.abs)),xl=sys=="ss"?a:L;
   g=`<path d="M${X(0)} 80H${X(L)}" stroke="#3b5f8a" stroke-dasharray="4 4"/><path d="M${xs.map((x,i)=>`${i?"L":""}${X(x)} ${80-Y[i]/dm*28}`).join("")}" fill="none" stroke="#38bdf8" stroke-width="3"/>`+(sys=="ss"?tri(X(0))+tri(X(L)):`<rect x="30" y="58" width="10" height="44" fill="#2a3f63"/>`)+arr(X(xl),22,72)+`<text x="${X(xl)+9}" y="34" fill="#fbbf24" font-size="12">${p} kN</text>`+
   `<text x="40" y="128" fill="#8da2bf" font-size="11">SHEAR FORCE (max ${f2(vm)} kN)</text>${path(xs,V,175,40/vm,"#38bdf8")}<text x="40" y="228" fill="#8da2bf" font-size="11">BENDING MOMENT, sagging down (max ${f2(mm)} kN·m)</text>${path(xs,M,270,-40/mm,"#34d399")}`;
   o=info.map(([k,v])=>`<div>${k}<b>${v.replace(/^-?\d+(\.\d+)?/,m=>f2(+m))}</b></div>`).join("")+`<div>δ_max<b>${f2(dm*1000)} mm</b></div>`}
  svg.innerHTML=g;out.innerHTML=o}
 [$s,$p,$a].forEach(e=>e.oninput=run);run()})();
/* ---- process animation ---- */
(function(){const el=qa("#proc [data-s]"),st=qa("#steps li"),bar=q("#pbar");let s=0,tm;const set=n=>{s=n;el.forEach(e=>e.style.opacity=+e.dataset.s<=n?1:0);st.forEach((e,i)=>e.classList.toggle("on",i==n-1));bar.style.width=n>=4?"380px":"0px"};
 const go=()=>set(s>=5?0:s+1);set(RM?5:0);if(!RM&&"IntersectionObserver"in window)new IntersectionObserver(e=>{clearInterval(tm);if(e[0].isIntersecting)tm=setInterval(go,2200)}).observe(q("#proc"))})();
/* ---- workspace preview: real numbers for a simply supported beam, UDL + midspan point load ---- */
(function(){const L=7,w=12,P=25,EI=40000,n=141,xs=[...Array(n).keys()].map(i=>i*L/(n-1)),R=(w*L+P)/2;
 const V=xs.map(x=>R-w*x-(x>L/2?P:0)),M=xs.map(x=>R*x-w*x*x/2-(x>L/2?P*(x-L/2):0)),Mx=Math.max(...M),dl=5*w*L**4/(384*EI)+P*L**3/(48*EI),X=x=>30+x/L*470;
 const pl=(v,y0,sc,c)=>`<path d="M${X(0)} ${y0} ${xs.map((x,i)=>`L${X(x)} ${y0-v[i]*sc}`).join(" ")} L${X(L)} ${y0}Z" fill="${c}33" stroke="${c}" stroke-width="1.5"/>`;
 q("#ws-svg").innerHTML=`<path d="M${X(0)} 70H${X(L)}" stroke="#8aa4d6" stroke-width="8"/>`+Array.from({length:15},(_,i)=>`<path d="M${X(i*L/14)} 28V62m-3 -7l3 7l3 -7" stroke="#fbbf24" fill="none"/>`).join("")+`<path d="M${X(L/2)} 8V62m-4 -8l4 8l4 -8" stroke="#fbbf24" stroke-width="3" fill="none"/><text x="${X(L/2)+8}" y="16" fill="#fbbf24" font-size="12">${P} kN</text><text x="${X(L)-110}" y="26" fill="#fbbf24" font-size="11">${w} kN/m</text><path d="M${X(0)} 74l-8 14h16zM${X(L)} 74l-8 14h16z" fill="none" stroke="#c9d6e8"/>`+
 `<text x="30" y="112" fill="#8da2bf" font-size="11">SHEAR FORCE DIAGRAM</text>${pl(V,165,45/Math.max(...V.map(Math.abs)),"#38bdf8")}<text x="30" y="226" fill="#8da2bf" font-size="11">BENDING MOMENT DIAGRAM (sagging down)</text>${pl(M,250,-45/Mx,"#34d399")}<text x="${X(L/2)+6}" y="308" fill="#34d399" font-size="11">M_max = ${f2(Mx)} kN·m</text>`;
 q("#ws-res").innerHTML=`<div>Reactions R_A = R_B<b>${f2(R)} kN</b></div><div>Max shear<b>${f2(Math.max(...V.map(Math.abs)))} kN</b></div><div>Max moment<b>${f2(Mx)} kN·m</b></div><div>Mid-span deflection<b>${f2(dl*1000)} mm</b></div>`})();
/* ---- reveal ---- */
const io=new IntersectionObserver(e=>e.forEach(x=>x.isIntersecting&&x.target.classList.add("in")),{threshold:.12});qa(".rv").forEach(e=>io.observe(e));
