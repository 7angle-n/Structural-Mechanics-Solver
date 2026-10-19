/* Determinate 2D pin-jointed truss. Tension +. Reactions solved by global equilibrium, members by method of joints. */
const gg=v=>(+v.toPrecision(4)).toString();
function gauss(A,b){const n=b.length;A=A.map((r,i)=>[...r,b[i]]);
 for(let c=0;c<n;c++){let p=c;for(let r=c+1;r<n;r++)if(Math.abs(A[r][c])>Math.abs(A[p][c]))p=r;if(Math.abs(A[p][c])<1e-9)throw Error("The truss is unstable or has a mechanism (singular system)");[A[c],A[p]]=[A[p],A[c]];
  for(let r=0;r<n;r++)if(r!=c){const k=A[r][c]/A[c][c];for(let j=c;j<=n;j++)A[r][j]-=k*A[c][j]}}
 return A.map((r,i)=>r[n]/r[i])}
function truss(){
 try{run()}catch(e){$("ans").innerHTML=`<p style="color:#c0392b">${e.message}</p>`;$("st").innerHTML="";$("svg").innerHTML=""}}
function run(){
 const rows=id=>$(id).value.trim().split("\n").map(l=>l.trim().split(/\s+/)).filter(a=>a[0]),idx={};
 const Nd=rows("nd").map((a,i)=>(idx[a[0]]=i,{n:a[0],x:+a[1],y:+a[2]}));
 const at=k=>{if(!(k in idx))throw Error("Unknown node "+k);return idx[k]};
 const M=rows("mb").map(a=>{const i=at(a[0]),j=at(a[1]),dx=Nd[j].x-Nd[i].x,dy=Nd[j].y-Nd[i].y,L=Math.hypot(dx,dy);if(!L)throw Error("Zero-length member "+a[0]+a[1]);return{i,j,L,c:dx/L,s:dy/L,nm:a[0]+a[1]}});
 const R=[];rows("sp").forEach(a=>{const i=at(a[0]);if(a[1]=="pin")R.push({i,d:0,nm:`R_{${a[0]}x}`},{i,d:1,nm:`R_{${a[0]}y}`});else{const d=a[2]=="x"?0:1;R.push({i,d,nm:`R_{${a[0]}${a[2]||"y"}}`})}});
 const P=Nd.map(()=>[0,0]);rows("ld").forEach(a=>{const i=at(a[0]);P[i][0]+=+a[1];P[i][1]+=+a[2]});
 const m=M.length,r=R.length,j=Nd.length;if(m+r!=2*j)throw Error(`Not statically determinate: m + r = ${m}+${r} = ${m+r}, but 2j = ${2*j}.`);
 const A=Array.from({length:2*j},()=>Array(m+r).fill(0)),b=[];
 M.forEach((q,k)=>{A[2*q.i][k]+=q.c;A[2*q.i+1][k]+=q.s;A[2*q.j][k]-=q.c;A[2*q.j+1][k]-=q.s});R.forEach((q,k)=>A[2*q.i+q.d][m+k]=1);
 P.forEach(p=>b.push(-p[0],-p[1]));const u=gauss(A,b),F=u.slice(0,m),Rv=u.slice(m);
 const gx=P.reduce((s,p)=>s+p[0],0),gy=P.reduce((s,p)=>s+p[1],0),O=Nd[R[0].i];
 const L=[{h:"Global equilibrium: reactions"},`\\sum F_x=0:\\ ${R.filter(q=>!q.d).map(q=>q.nm).join("+")||"0"}+(${gx})=0`,`\\sum F_y=0:\\ ${R.filter(q=>q.d).map(q=>q.nm).join("+")||"0"}+(${gy})=0`,
  `\\sum M_{${O.n}}=0:\\ ${R.map(q=>`(${gg(q.d?Nd[q.i].x-O.x:-(Nd[q.i].y-O.y))})${q.nm}`).join("+")}+(${gg(P.reduce((s,p,i)=>s+(Nd[i].x-O.x)*p[1]-(Nd[i].y-O.y)*p[0],0))})=0`,
  R.map((q,k)=>`${q.nm}=${gg(Rv[k])}`).join(",\\ ")];
 const done=Array(m).fill(false),jd=Array(j).fill(false);let prog=true;L.push({h:"Method of joints (tension +, compression −)"});
 while(prog){prog=false;for(let n=0;n<j;n++){if(jd[n])continue;
  const un=M.map((q,k)=>[q,k]).filter(([q,k])=>(q.i==n||q.j==n)&&!done[k]);if(!un.length){jd[n]=true;continue}if(un.length>2)continue;
  const co=([q])=>q.i==n?[q.c,q.s]:[-q.c,-q.s],K=[P[n][0],P[n][1]];
  R.forEach((q,k)=>{if(q.i==n)K[q.d]+=Rv[k]});M.forEach((q,k)=>{if(done[k]&&(q.i==n||q.j==n)){const c=co([q]);K[0]+=c[0]*F[k];K[1]+=c[1]*F[k]}});
  let sol;if(un.length==2){const[a,c]=[co(un[0]),co(un[1])],det=a[0]*c[1]-c[0]*a[1];if(Math.abs(det)<1e-9)continue;sol=[(-K[0]*c[1]+K[1]*c[0])/det,(-a[0]*K[1]+a[1]*K[0])/det]}
  else{const a=co(un[0]),ax=Math.abs(a[0])>Math.abs(a[1])?0:1;sol=[-K[ax]/a[ax]]}
  const ex=[0,1].map(ax=>un.map(w=>`(${gg(co(w)[ax])})F_{${w[0].nm}}`).join("+")+`+(${gg(K[ax])})=0`);
  L.push({h:"Joint "+Nd[n].n},`\\sum F_x=0:\\ ${ex[0]}`,`\\sum F_y=0:\\ ${ex[1]}`,un.map((w,i)=>`F_{${w[0].nm}}=${gg(sol[i])}`).join(",\\ "));
  un.forEach((w,i)=>{done[w[1]]=true;F[w[1]]=sol[i]});jd[n]=true;prog=true}}
 if(done.includes(false))L.push({h:"Remaining members"},"\\text{Solved simultaneously from the full equilibrium system (joint-by-joint order stalls).}");
 steps($("st"),L);
 const Tm=Math.max(0,...F),Cm=Math.min(0,...F);
 $("ans").innerHTML=R.map((q,k)=>`<div>${q.nm.replace(/_\{(.*)\}/,"$1")}<b>${gg(Rv[k])}</b></div>`).join("")+`<div>Max tension<b>${gg(Tm)}</b></div><div>Max compression<b>${gg(Cm)}</b></div>`;
 $("tb").innerHTML="<tr><th>Member</th><th>Length</th><th>Force</th><th></th></tr>"+M.map((q,k)=>`<tr><td>${q.nm}</td><td>${gg(q.L)}</td><td>${gg(F[k])}</td><td>${Math.abs(F[k])<1e-9?"zero-force":F[k]>0?"Tension":"Compression"}</td></tr>`).join("");
 const xs=Nd.map(p=>p.x),ys=Nd.map(p=>p.y),x0=Math.min(...xs),y0=Math.min(...ys),k=Math.min(480/((Math.max(...xs)-x0)||1),220/((Math.max(...ys)-y0)||1)),X=v=>50+(v-x0)*k,Y=v=>270-(v-y0)*k;
 let g="";M.forEach((q,i)=>{const a=Nd[q.i],c=Nd[q.j],col=Math.abs(F[i])<1e-9?"#999":F[i]>0?"#1d4ed8":"#c0392b";g+=`<line x1="${X(a.x)}" y1="${Y(a.y)}" x2="${X(c.x)}" y2="${Y(c.y)}" stroke="${col}" stroke-width="4" stroke-linecap="round"/><text x="${(X(a.x)+X(c.x))/2+6}" y="${(Y(a.y)+Y(c.y))/2-6}" fill="${col}" font-size="15">${gg(F[i])}</text>`});
 R.forEach(q=>{const p=Nd[q.i];g+=`<path d="M${X(p.x)} ${Y(p.y)}l-10 18h20z" fill="none" stroke="#33413d" stroke-width="2"/>`});
 P.forEach((p,i)=>{if(p[0]||p[1]){const n=Math.hypot(...p),dx=p[0]/n*36,dy=-p[1]/n*36;g+=`<path d="M${X(Nd[i].x)-dx} ${Y(Nd[i].y)-dy}L${X(Nd[i].x)} ${Y(Nd[i].y)}" stroke="#a0612b" stroke-width="3" marker-end="url(#ah)"/>`}});
 Nd.forEach(p=>g+=`<circle cx="${X(p.x)}" cy="${Y(p.y)}" r="5" fill="#33413d"/><text x="${X(p.x)-14}" y="${Y(p.y)+22}" font-size="16" fill="#33413d">${p.n}</text>`);
 window.TR={Nd,M,R,Rv,P,F};$("svg").innerHTML=`<defs><marker id="ah" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0 0L8 4L0 8z" fill="#a0612b"/></marker></defs>`+g;
}

function sections(){try{sect()}catch(e){$("ans2").innerHTML=`<p style="color:#c0392b">${e.message}</p>`;$("st2").innerHTML=""}}
function sect(){
 if(!window.TR)truss();const{Nd,M,R,Rv,P,F}=window.TR,ks=$("cut").value.trim().split(/\s+/).map(n=>{const k=M.findIndex(q=>q.nm==n||q.nm==[...n].reverse().join(""));if(k<0)throw Error("Unknown member "+n);return k});
 if(ks.length!=3)throw Error("Choose exactly 3 members to cut");
 const adj=Nd.map(()=>[]);M.forEach((q,k)=>{if(!ks.includes(k)){adj[q.i].push(q.j);adj[q.j].push(q.i)}});
 let s0=0;Nd.forEach((p,i)=>{if(p.x<Nd[s0].x||(p.x==Nd[s0].x&&p.y<Nd[s0].y))s0=i});
 const inn=new Set([s0]),st=[s0];while(st.length){const u=st.pop();adj[u].forEach(v=>{if(!inn.has(v)){inn.add(v);st.push(v)}})}
 const cut=ks.map(k=>{const q=M[k],ins=inn.has(q.i);if(ins==inn.has(q.j))throw Error("Those members do not separate the truss into two parts");const n=ins?q.i:q.j,o=ins?q.j:q.i,dx=Nd[o].x-Nd[n].x,dy=Nd[o].y-Nd[n].y,l=Math.hypot(dx,dy);return{k,n,c:dx/l,s:dy/l,nm:q.nm}});
 const O=Nd[s0];let Kx=0,Ky=0,Km=0;inn.forEach(n=>{Kx+=P[n][0];Ky+=P[n][1];Km+=(Nd[n].x-O.x)*P[n][1]-(Nd[n].y-O.y)*P[n][0]});
 R.forEach((q,k)=>{if(inn.has(q.i)){const f=Rv[k];if(q.d==0){Kx+=f;Km+=-(Nd[q.i].y-O.y)*f}else{Ky+=f;Km+=(Nd[q.i].x-O.x)*f}}});
 const mo=cut.map(q=>(Nd[q.n].x-O.x)*q.s-(Nd[q.n].y-O.y)*q.c);let f;
 try{f=gauss([cut.map(q=>q.c),cut.map(q=>q.s),mo],[-Kx,-Ky,-Km])}catch(e){throw Error("The three cut members are concurrent or parallel; choose another section")}
 const row=(co,K)=>cut.map((q,i)=>`(${gg(co[i])})F_{${q.nm}}`).join("+")+`+(${gg(K)})=0`;
 steps($("st2"),[{h:"Free body: nodes "+[...inn].map(i=>Nd[i].n).join(", ")},`\\text{Cut members }${cut.map(q=>q.nm).join(",\\ ")}\\text{; assume tension (force pulls the free body toward the cut).}`,
  `\\sum F_x=0:\\ ${row(cut.map(q=>q.c),Kx)}`,`\\sum F_y=0:\\ ${row(cut.map(q=>q.s),Ky)}`,`\\sum M_{${O.n}}=0:\\ ${row(mo,Km)}`,
  cut.map((q,i)=>`F_{${q.nm}}=${gg(f[i])}`).join(",\\ "),`\\text{Check against the joint method: }${cut.map((q,i)=>`${q.nm}:\\ ${gg(F[q.k])}`).join(",\\ ")}`]);
 $("ans2").innerHTML="";$("ans2").innerHTML=cut.map((q,i)=>`<div style="display:inline-block;margin-right:22px">${q.nm}<b style="font-size:1.6rem">${gg(f[i])}</b>${f[i]>1e-9?"tension":f[i]<-1e-9?"compression":"zero"}</div>`).join("")}
